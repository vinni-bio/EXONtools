# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import print_function, division
import os
import logging
import csv
import pdb

from mains.EXT_prog import EXTprogram
from mains.EXT_IO import getinput, output, makenewdir
from mains.EXT_executor import executor
from mains.EXT_parallel import run_executor
from mains.EXT_errors import EXONtoolsError
from progs.annotator import test_format
from utils.revcomp import DNArevcomp
from utils.sorting import natural_sort
from utils.seqIO import SeqIO
from utils.blastout import BLAST6


class exosearcher(EXTprogram):
    """This program annotates the provided assembly using the annotation reference"""

    name = "exosearcher"

    def execute_program(self):
        args = self.args
        self.search_exons(args.inpath, args.gff, args.reference, args.outdir, args.evalue, args.lag, args.similarity, args.overlap, args.chimeric, args.isomeric, args.naexclude, args.unique, args.orfcheck)

    def search_exons(self, inpath, gffpath, reference, outdir, evalue, lag, similarity, overlap, chimeric, isomeric, naexclude, unique, orfcheck):

        if exosearcher.debug:
            pdb.set_trace()

        executor.setconfig("makeblastdb", "blastn")

        if exosearcher.debug:
            pdb.set_trace()

        exosearcher_pars(evalue, lag, similarity, overlap, chimeric, isomeric, naexclude, unique, orfcheck)

        # GET INPUT FILES
        # Reference file checkup
        getinput.format(['.fa', '.fasta', '.fna'])
        reffiles = getinput(reference).files
        FileList = getinput(inpath).files
        getinput.format(['.gff', '.gff3', '.gtf'])
        AnnotList = getinput(gffpath).files

        # MAKE OUTPUT DIRECTORY
        output(outdir)

        # MAKE TMP DIRECTORIES
        tmpdir = makenewdir(name="tmp", fullname="temporary")

        # MAKE LOG DIR
        log_dir = makenewdir(name="dependency_logs", fullname="LOG")
        logging.info("All dependency program logs will be saved in the '{0:s}' folder".format(os.path.basename(log_dir.path)))

        logging.debug("IO settings: OK")

        logging.debug("Ready to start the annotation analysis")

        if exosearcher.debug:
            pdb.set_trace()

        stats_collector = {}

        logging.info("**********************************************************************")
        logging.info("STAGE 1. PREPARING ASSEMBLY ANNOTATIONS")

        FileList.sort(key=natural_sort)
        AnnotList.sort(key=natural_sort)
        LibDict = {}
        for x, y in zip(FileList, AnnotList):
            xfile = os.path.basename(x)
            yfile = os.path.basename(y)
            lib = xfile.split("_")[0].split(".")[0]
            if lib != yfile.split("_")[0].split(".")[0]:
                logging.error("Assembly and annotation files do not match by their library identifiers!")
                raise EXONtoolsError("Input file format error")
            else:
                LibDict[lib] = (x, y)
        logging.debug("All input assembly files have corresponding annotation files: OK")
        logging.info("Total {0:d} annotated assemblies are provided for exon search analysis".format(len(LibDict)))

        if exosearcher.debug:
            pdb.set_trace()

        contigs = {}
        selected = makenewdir(name=os.path.join(tmpdir.path, "SELECTED"), fullname="SELECTED")

        for lib in sorted(LibDict.keys(), key=natural_sort):
            logging.info("Parsing annotation in '{0:s}' library".format(lib))

            y = test_format([LibDict[lib][0]], stype="nucl")
            stats_collector[lib] = {"Total": y[lib]}
            selected_annots = {}
            contigs[lib] = {}
            filtered = []
            chimercount = 0
            selpaths = {}

            annotoutpath = os.path.join(output.path, lib + "_reference.gff")

            with open(LibDict[lib][1], 'r') as gffile:
                for line in gffile:
                    line = line.strip()
                    if line.startswith("##sequence-region"):
                        llist = line.split("\t")
                        contigs[lib][llist[1].strip()] = {"length": int(llist[3].strip())}
                    elif not line or line.startswith("#"):
                        pass
                    else:
                        llist = line.split("\t")
                        if llist[2].strip() == 'gene':
                            contigname = llist[0].strip()
                            contigs[lib][contigname]['name'] = contigname
                            infolist = llist[8].split(";")
                            isoforms = []
                            ischimeric = False
                            for info in infolist:
                                info = info.strip()
                                if info.startswith("ID="):
                                    geneid = info.replace("ID=", "")
                                elif info.startswith("isoforms=") and isomeric:
                                    isoforms = info.replace("isoforms=", "").split(",")
                                elif info.startswith("Name="):
                                    genename = info.replace("Name=", "")
                                elif info.startswith("is_chimeric=true") and chimeric:
                                    ischimeric = True
                                elif info.startswith("gene="):
                                    gene = info.replace("gene=", "").strip().lower()
                                    if gene.upper() == 'NA' and (naexclude or unique):
                                        filtered.append(contigname)
                                        continue
                                    elif unique:
                                        try:
                                            if contigs[lib][contigname]['length'] >= contigs[lib][selected_annots[gene]]['length']:
                                                filtered.append(selected_annots[gene])
                                                # contigs[lib][selected_annots[gene]] = {contigs[lib][selected_annots[gene]]['length']}
                                                selected_annots[gene] = contigname
                                            else:
                                                filtered.append(contigname)
                                                continue
                                        except KeyError:
                                            selected_annots[gene] = contigname
                                    else:
                                        pass
                                else:
                                    pass

                            if isomeric and [x for x in isoforms if contigs[lib][x]["length"] >= contigs[lib][contigname]["length"]]:
                                filtered.append(contigname)
                            else:
                                contigs[lib][contigname]["annotation"] = genename
                                contigs[lib][contigname]["gene"] = gene
                                contigs[lib][contigname]["parent"] = geneid
                                contigs[lib][contigname]["chimeric"] = ischimeric
                                contigs[lib][contigname]["start"] = int(llist[3])
                                contigs[lib][contigname]["end"] = int(llist[4])
                                contigs[lib][contigname]["strand"] = llist[6]

                        if llist[2].strip() == 'CDS':
                            contigname = llist[0].strip()
                            start = int(llist[3].strip())
                            end = int(llist[4].strip())
                            strand = llist[6].strip()
                            frame = llist[7].strip()
                            contigs[lib][contigname]['CDS'] = (start, end, strand, frame)

            filtered = list(set(filtered))
            for x in filtered:
                del contigs[lib][x]

            if not contigs[lib]:
                logging.error("All contigs in '{0:s}' library were filtered. Please check command parsing settings".format(lib))
                raise EXONtoolsError("No contigs were selected for the exon prediction analysis")

            outpath = os.path.join(selected.path, lib + "_selected.fasta")
            selpaths[lib] = outpath
            with open(outpath, 'w') as outfile:
                fastafile = SeqIO(LibDict[lib][0])
                for contig in fastafile.read():
                    try:
                        name = contigs[lib][contig.name.split("\t")[0]]['name']
                        seq = contig.seq
                        if contigs[lib][name]["chimeric"]:
                            chimercount += 1
                            scut = max(0, contigs[lib][name]["start"] - lag)
                            ecut = min(contigs[lib][name]["length"], contigs[lib][name]["end"] + lag + 1)
                            seq = seq[scut:ecut]
                        if contigs[lib][name]["strand"] == "-":
                            seq = DNArevcomp(seq)
                        outfile.write(">{0:s}\n{1:s}\n".format(name, seq))
                    except KeyError:
                        pass

            stats_collector[lib]['filtered'] = len(filtered)
            stats_collector[lib]['chimeric'] = chimercount
            logging.info("{0:d} contigs were filtered and {1:d} chimeras were corrected in '{2:s}' library".format(len(filtered), chimercount, lib))

        if exosearcher.debug:
            pdb.set_trace()

        logging.info("**********************************************************************")
        logging.info("STAGE 2. PREPARING REFERENCE GENOMES")
        logging.info("Total {0:d} reference genomes are provided for exon search analysis".format(len(reffiles)))
        REFDIR = makenewdir(name=os.path.join(tmpdir.path, "REFERENCE"), fullname="REFERENCE")
        refpath = os.path.join(REFDIR.path, "REFS.fasta")
        refdbpath = os.path.join(REFDIR.path, "REFS")
        test_format(reffiles, stype="nucl")
        logging.debug("Writing all references into one file to construct BLAST database")
        with open(refpath, 'w') as outfile:
            for ref in reffiles:
                with open(ref, 'r') as infile:
                    for line in infile:
                        outfile.write(line)
        outlog = os.path.join(log_dir.path, "blast.log")
        logging.info("Constructing reference database")
        run_executor(executor(
            program="makeblastdb",
            params=[refpath, refdbpath, 'nucl'],
            conditions={"pathexists": [refpath]},
            custom_arg_string=" >> " + outlog + " 2>&1",
            quiet=True)
        )

        if not exosearcher.dryrun:
            os.remove(refpath)
        logging.debug("BLAST database is constructed: OK")

        if exosearcher.debug:
            pdb.set_trace()

        logging.info("**********************************************************************")
        logging.info("STAGE 3. MAPPING CONTIGS TO REFERENCES")
        BLASTDIR = makenewdir(name=os.path.join(tmpdir.path, "BLAST"), fullname="BLAST")
        blastpars = "-word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -template_type coding -template_length 18"
        blastpaths = {}

        for lib in sorted(selpaths.keys(), key=natural_sort):
            blastpath = os.path.join(BLASTDIR.path, lib + ".outfmt6")
            blastpaths[lib] = blastpath
            logging.info("Mapping contigs in '{0:s}' library".format(lib))

            run_executor(executor(
                program="blastn",
                params=[refdbpath, selpaths[lib], blastpath, 6, str(evalue), exosearcher.threads, blastpars],
                conditions={"pathexists": [selpaths[lib]], "positive": [exosearcher.threads]},
                custom_arg_string=" >> " + outlog + " 2>&1")
            )
        logging.debug("BLAST alignment is finished: OK")

        if exosearcher.debug:
            pdb.set_trace()

        logging.info("**********************************************************************")
        logging.info("STAGE 4. PREDICTING EXON BOUNDARIES")

        logging.debug("Collecting BLAST results")

        BLASTRES = {}

        for lib in sorted(blastpaths.keys(), key=natural_sort):
            logging.info("Predicting exon boundaries for '{0:s}' library".format(lib))
            BLASTRES[lib] = {}
            blastout = BLAST6(blastpaths[lib])
            previous = None
            for bline in blastout.read():
                if bline.length >= overlap:
                    if bline.tend - bline.tstart >= 0:
                        sign = "+"
                    else:
                        sign = "-"
                    try:
                        BLASTRES[lib][bline.query].append((bline.query, bline.target, bline.identity, sign, bline.qstart, bline.qend, bline.evalue, bline.bit, bline.tstart, bline.tend))
                    except KeyError:
                        BLASTRES[lib][bline.query] = [(bline.query, bline.target, bline.identity, sign, bline.qstart, bline.qend, bline.evalue, bline.bit, bline.tstart, bline.tend)]
                        try:
                            BLASTRES[lib][previous] = finalize_exons(previous, BLASTRES[lib][previous],contigs[lib][previous],lag,similarity)
                        except KeyError:
                            pass
                        previous = bline.query
            try:
                BLASTRES[lib][previous] = finalize_exons(previous, BLASTRES[lib][previous],contigs[lib][previous],lag,similarity)
            except KeyError:
                pass

            #### CLEANING RESULTS
            for contig in list(BLASTRES[lib].keys()):
                if isinstance(BLASTRES[lib][contig],list):
                    logging.error("Not all blast output was parsed correctly")
                    raise EXONtoolsError("BLAST parsing error")
                elif not BLASTRES[lib][contig]:
                    del BLASTRES[lib][contig]
                else:
                    pass

            stats_collector[lib]['blast'] = len(BLASTRES[lib])

        logging.debug("Exon prediction sucessfully completed: OK")

        if exosearcher.debug:
            pdb.set_trace()

        logging.info("**********************************************************************")
        logging.info("STAGE 5. SAVING EXON PREDICTION RESULTS")

        exon_lengths = {}
        for lib in sorted(BLASTRES.keys(), key=natural_sort):
            exon_lengths[lib] = []
            outexonpath = os.path.join(output.path, lib + exosearcher.suffix + ".fasta")
            outfastapath = os.path.join(output.path, lib + "_reference.fasta")

            logging.info("Writing exons and reference files for '{0:s}' library".format(lib))
            infastafile = SeqIO(LibDict[lib][0])
            outexonfile = open(outexonpath, 'w')
            outfastafile = open(outfastapath, 'w')

            for contig in infastafile.read():
                contigname = contig.name.split("\t")[0]
                try:
                    predictions = BLASTRES[lib][contigname]
                    outfastafile.write(">{0:s}\n{1:s}\n".format(contig.name, contig.seq))
                except KeyError:
                    continue

                if orfcheck and 'CDS' in contigs[lib][contigname]:
                    checkrange = set(range(contigs[lib][contigname]['CDS'][0], contigs[lib][contigname]['CDS'][1] + 1))
                else:
                    checkrange = None

                for exon in sorted(predictions.keys(), key=natural_sort):
                    exoseq = contig.seq[predictions[exon][4] - 1:predictions[exon][5]]
                    exon_lengths[lib].append(len(exoseq))
                    if checkrange:
                        if checkrange & set(range(predictions[exon][4], predictions[exon][5] + 1)):
                            outexonfile.write(">{0:s}_{1:s}_{2:d}_{3:d}\n{4:s}\n".format(contigname, exon, predictions[exon][4], predictions[exon][5], exoseq))
                    else:
                        outexonfile.write(">{0:s}_{1:s}_{2:d}_{3:d}\n{4:s}\n".format(contigname, exon, predictions[exon][4], predictions[exon][5], exoseq))

            outexonfile.close()
            outfastafile.close()

            annotoutfile = open(annotoutpath, 'w')
            with open(LibDict[lib][1], 'r') as gffile:
                exoncount = 0
                for line in gffile:
                    line = line.strip()
                    if line.startswith("##sequence-region"):
                        llist = line.split("\t")
                        if llist[1].strip() in BLASTRES[lib]:
                            annotoutfile.write(line + "\n")
                    elif line.startswith("#"):
                        annotoutfile.write(line + "\n")
                    else:
                        llist = line.split("\t")
                        contigname = llist[0].strip()
                        if contigname in BLASTRES[lib]:
                            if llist[2].strip() == 'gene':
                                annotoutfile.write(line + "\n")
                                predictions = BLASTRES[lib][contigname]
                                for exon in sorted(predictions.keys(), key=natural_sort):
                                    exoncount += 1
                                    exonatribs = "ID=exon{0:d};Parent={1:s};Target={2:s} {3:d} {4:d} {5:s}".format(exoncount, contigs[lib][contigname]["parent"], predictions[exon][1], predictions[exon][8], predictions[exon][9], predictions[exon][3])
                                    annotoutfile.write("{0:s}\t{1:s}\texon\t{2:d}\t{3:d}\t{4:.2e}\t{5:s}\t.\t{6:s}\n".format(contigname, "EXONtools", predictions[exon][4], predictions[exon][5], predictions[exon][6], llist[6], exonatribs))
                            else:
                                annotoutfile.write(line + "\n")
                        else:
                            continue
            annotoutfile.close()

        logging.debug("All exons were successfully saved to files: OK")

        if exosearcher.debug:
            pdb.set_trace()

        if exosearcher.stats:
            logging.info("**********************************************************************")
            logging.info("STAGE 6. SAVING STATS FOR EXON PREDICTION ANALYSIS")
            logging.info("All exon prediction stats will be saved to 'exon_prediction_stats.csv' file")
            logging.info("Writing exon number distribution to 'lib_exon_numbs.txt'")
            if not exosearcher.dryrun:
                for lib in BLASTRES:
                    outpath = os.path.join(output.path, lib + "_exon_numbs.txt")
                    exon_nums = list(map(lambda x: len(BLASTRES[lib][x]),BLASTRES[lib]))
                    exon_nums.sort()
                    with open(outpath, 'w') as outfile:
                        for x in exon_nums:
                            outfile.write("{0:d}\n".format(x))
                logging.debug("Exon numbers per contig were successfully saved to files: OK")
            logging.info("Writing exon length distribution to 'lib_exon_lengths.txt'")
            if not exosearcher.dryrun:
                for lib in exon_lengths:
                    outpath = os.path.join(output.path, lib + "_exon_lengths.txt")
                    with open(outpath, 'w') as outfile:
                        for exval in sorted(exon_lengths[lib]):
                            outfile.write("{0:d}\n".format(exval))
                logging.debug("Exon lengths were successfully saved to files: OK")
            if not exosearcher.dryrun:
                logging.info("Writing stats...")
                statpath = os.path.join(output.path, "exon_prediction_stats.csv")
                header = ["No", "RRL", "N_CONTIGS", "FILTERED", "CHIMERAS", "MAPPED", "EXONS", "EXONS>100bp", "EXONS>200bp", "EXONS>500bp", "EXONS>1000bp", "EXONS>5000bp"]
                with open(statpath, 'w') as statfile:
                    csv_writer = csv.writer(statfile)
                    csv_writer.writerow(header)
                    for i, lib in enumerate(sorted(stats_collector.keys(), key=natural_sort)):
                        statline = [i + 1, lib, stats_collector[lib]['Total'], stats_collector[lib]['filtered'], stats_collector[lib]['chimeric'], stats_collector[lib]['blast']]
                        check100 = 0
                        check200 = 0
                        check500 = 0
                        check1000 = 0
                        check5000 = 0
                        for x in sorted(exon_lengths[lib], reverse=True):
                            if x > 5000:
                                check5000 += 1
                                continue
                            if x > 1000:
                                check1000 += 1
                                continue
                            if x > 500:
                                check500 += 1
                                continue
                            if x > 200:
                                check200 += 1
                                continue
                            if x > 100:
                                check100 += 1
                                continue
                        statline.append(len(exon_lengths[lib]))
                        statline.append(check100 + check200 + check500 + check1000 + check5000)
                        statline.append(check200 + check500 + check1000 + check5000)
                        statline.append(check500 + check1000 + check5000)
                        statline.append(check1000 + check5000)
                        statline.append(check5000)
                        csv_writer.writerow(statline)

        if not exosearcher.keeptmp:
            tmpdir.delete()


def finalize_exons(contig, blastlist,contdic,lag,similarity):
    """Predicting and fixing exon boundaries from blast output"""

    temp = {}
    EXONS = {}
    check = ("", 0, 0, 0)
    a = 0

    #### PARSING BLAST OUTPUT ####
    for i, x in enumerate(sorted(blastlist, key=lambda e: e[4])):

        temp["match" + str(i)] = x

        if x[4] > check[1] + lag and x[5] > check[2] + lag and x[2] > similarity:
            if check != ("", 0, 0, 0):
                a += 1
                EXONS["Exon" + str(a)] = temp[check[0]]
            check = ("match" + str(i), x[4], x[5], x[7])
        elif x[4] < lag and check == ("", 0, 0, 0):
            check = ("match" + str(i), x[4], x[5], x[7])
        else:
            if x[2] > similarity and x[7] > check[3]:
                check = ("match" + str(i), x[4], x[5], x[7])
            elif x[6] == 0 and x[7] > check[3]:
                check = ("match" + str(i), x[4], x[5], x[7])
            else:
                continue

    if check[0] != "":
        a += 1
        EXONS["Exon" + str(a)] = temp[check[0]]
    else:
        return None

    #### FILTERING AND CORRECTING STRAND DIRECTION####
    temps = []
    for exon in sorted(EXONS.keys(), key=natural_sort):
        exomatch = list(EXONS[exon])
        if contdic["strand"] == "-" and not contdic["chimeric"]:
            check = exomatch[4]
            exomatch[4] = contdic['length'] - exomatch[5] + 1
            exomatch[5] = contdic['length'] - check + 1
        elif contdic["strand"] == "+" and contdic["chimeric"]:
            scut = max(0, contdic["start"] - lag)
            exomatch[4] = exomatch[4] + scut
            exomatch[5] = exomatch[5] + scut
        elif contdic["strand"] == "-" and contdic["chimeric"]:
            check = exomatch[4]
            ecut = min(contdic["length"], contdic["end"] + lag + 1)
            exomatch[4] = ecut - exomatch[5] + 1
            exomatch[5] = ecut - check + 1
        else:
            pass
        temps.append(tuple(exomatch))

    EXONS = {}
    excounter = 0
    for temp in sorted(temps, key=lambda x: x[4]):
        if temp[5] - temp[4] >= lag:
            excounter += 1
            EXONS["Exon" + str(excounter)] = temp
    if not EXONS:
        return None

    #### CORRECT OVERLAPS ####
    temps = []
    previous = None
    for exon in sorted(EXONS.keys(), key=natural_sort):
        exomatch = list(EXONS[exon])
        if previous:

            if len(set(list(range(exomatch[4], exomatch[5] + 1))) - set(list(range(previous[4], previous[5] + 1)))) < 10:
                continue
            elif len(set(list(range(previous[4], previous[5] + 1))) - set(list(range(exomatch[4], exomatch[5] + 1)))) < 10:
                previous = list(temps[-1])[::]
                temps = temps[:-1]
            else:
                pass

            if previous[5] >= exomatch[4] and previous[6] < exomatch[6]:
                exomatch[4] = previous[5] + 1
            elif previous[5] >= exomatch[4] and previous[6] >= exomatch[6]:
                previous[5] = exomatch[4] - 1
            elif previous[5] < exomatch[4] and exomatch[4] - previous[5] < 10:
                if previous[6] < exomatch[6]:
                    exomatch[4] = previous[5] + 1
                else:
                    previous[5] = exomatch[4] - 1
            elif previous[5] < exomatch[4] and lag<exomatch[4] - previous[5]<200:
                temps.append(tuple([contig,"NA",0,exomatch[3], previous[5]+1,exomatch[4]-1,1,0,0,0]))
            else:
                pass
            temps.append(tuple(previous))
        previous = exomatch[::]

    temps.append(tuple(previous))
    EXONS = {}
    for i, temp in enumerate(sorted(temps, key=lambda x: x[4])):
        if temp[4]<temp[5]:
            EXONS["Exon" + str(i + 1)] = temp

    return EXONS


def exosearcher_pars(evalue, minlen, similarity, overlap, chimeric, isomeric, naexclude, unique, orfcheck):
    """Testing exosearcher parameters"""

    logging.debug("Testing exosearcher parameters")

    if evalue < 0:
        logging.error("The e-value argument cannot be negative")
        raise EXONtoolsError("Argument value error")
    else:
        logging.info("E-value alignment threshold is set to '{0:s}'".format(str(evalue)))

    if minlen < 0:
        logging.error("Minumum contig length cannot have a negative value")
        raise EXONtoolsError("Argument value error")

    if similarity < 0 or similarity > 1:
        logging.error("Alignment similarity value cannot be lower than zero or greater than 1 ")
        raise EXONtoolsError("Argument value error")

    if overlap < 0:
        logging.error("Alignment length threshold cannot be less than 0")
        raise EXONtoolsError("Argument value error")

    if chimeric:
        logging.warning("Chimeric contigs will be parsed to choose the correctly annotated region for mapping")

    if isomeric:
        logging.warning("Contig isomers will be filtered. Only the longest isomer will be used for exon prediction")
    else:
        logging.warning("All contig isomers will be analyzed for exon prediction")

    if unique:
        logging.warning("Only the longest contig will be used for each annotation (unique annotations in the output)")

    if naexclude or unique:
        logging.warning("All contigs with unknown annotation (e.g. gene=NA) will be excluded from the analysis")

    if orfcheck:
        logging.warning("Only exons that overlap with the CDS regions will be provided in the final output fasta file")
