# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

from __future__ import print_function, division
import os
import re
import logging
import csv
import pdb
from decimal import Decimal

from mains.EXT_prog import EXTprogram
from mains.EXT_IO import getinput, output, makenewdir
from mains.EXT_executor import executor
from mains.EXT_worker import worker
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, run_executor, set_threads
from mains.EXT_errors import EXONtoolsError
from utils.sorting import natural_sort
from utils.seqIO import SeqIO
from utils.tandreps import microsat
from utils.encoders import findorfs, translate, loadcodes
from utils.blastout import BLAST6
from utils.revcomp import DNArevcomp


class annotator(EXTprogram):
    """This program annotates the provided assembly using the annotation reference"""

    name = "annotator"

    def execute_program(self):
        args = self.args
        self.annotate_seqs(args.inpath, args.reference, args.filtering, args.scaffolds, args.isoforms, args.chimeras, args.outdir, args.query, args.target, args.evalue, args.cluster, args.minlen, args.similarity, args.overlap, args.database, args.custom, args.orf, args.translating, args.gencode,args.oneway, args.norename, args.organism)

    def annotate_seqs(self, inpath, reference, filtering, scaffolds, isoforms, chimeras, outdir, query, target, evalue, cluster, minlen, similarity, overlap, database, custom, orf, translating, gencode, oneway, norename,organism):

        if annotator.debug:
            pdb.set_trace()

        # TESTING INPUT PARAMETER VALUES
        mem = int(annotator.memory * 1000)

        # DICTIONARY FOR COLLECTING STATS
        annostats = {"targets": {}, "queries": {}}

        if annotator.debug:
            pdb.set_trace()
        executor.setconfig("makeblastdb", "blastn", "blastx", "blastp", "tblastx", "tblastn", "cd-hit-est")
        if annotator.debug:
            pdb.set_trace()
        annotator_pars(query, target, evalue, cluster, minlen, similarity, overlap, database, custom, orf, translating, gencode, oneway, norename, organism)
        if annotator.debug:
            pdb.set_trace()

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        annotator.run_dry(BLAST6, SeqIO, getinput, output, makenewdir, executor, annotat)
        annotator.set_debug(executor)

        # GET INPUT FILES
        # Reference file checkup
        getinput.format(['.fa', '.fasta', '.prot', '.pep'])
        if os.path.isfile(reference):
            refpath = getinput(reference).files[0]
        else:
            logging.error("Reference file not found. Please verify provided path and try again")
            raise EXONtoolsError("Input path error")
        FileList = getinput(inpath).files

        # GET SUPPLEMENTARY FILES
        getinput.setdry()
        if scaffolds:
            getinput.format(['.scaffolds'])
            if os.path.exists(scaffolds):
                scaffolds = {os.path.basename(x).split(".")[0]: x for x in getinput(scaffolds).files}
                if scaffolds:
                    logging.warning("Scaffolds are provided for ANNOTATION analysis")
                else:
                    logging.warning("No scaffolds were found in the provided file. Continue analysis...")
            else:
                logging.error("No scaffolds are found in the provided path")
                raise EXONtoolsError("Input path error")
        if isoforms:
            getinput.format(['.isoforms'])
            if os.path.exists(isoforms):
                isoforms = {os.path.basename(x).split(".")[0]: x for x in getinput(isoforms).files}
                if isoforms:
                    logging.warning("Isoforms are provided for ANNOTATION analysis")
                else:
                    logging.warning("No isoforms were found in provided file. Continue analysis...")
            else:
                logging.error("No isoforms are found in the provided path")
                raise EXONtoolsError("Input path error")
        if chimeras:
            getinput.format(['.chimeras'])
            if os.path.exists(chimeras):
                chimeras = {os.path.basename(x).split(".")[0]: x for x in getinput(chimeras).files}
                if chimeras:
                    logging.warning("Chimeras are provided for ANNOTATION analysis")
                else:
                    logging.warning("No chimeras were found in the provided file. Continue analysis...")
            else:
                logging.error("No chimeras are found in provided path")
                raise EXONtoolsError("Input path error")

        if not annotator.dryrun:
            getinput.undry()

        # MAKE OUTPUT DIRECTORY
        output(outdir)

        # MAKE TMP DIRECTORIES
        tmpdir = makenewdir(name="tmp", fullname="temporary")

        # MAKE LOG DIR
        log_dir = makenewdir(name="dependency_logs", fullname="LOG")
        logging.info("All dependency program logs will be saved in the '{0:s}' folder".format(os.path.basename(log_dir.path)))

        logging.debug("IO settings: OK")

        logging.debug("Ready to start the annotation analysis")

        if annotator.debug:
            pdb.set_trace()

        logging.info("**********************************************************************")
        logging.info("STAGE 1. REFERENCE DATABASE CONSTRUCTION")
        targresult = test_format([refpath], target)
        annostats["targets"]["total"] = [x for x in targresult.values()][0]
        del targresult

        # DEFINE TYPE OF REFERENCE FOR PARSING
        if database.lower() == "swissprot":
            ref_type = 'SWISSPROT'
        elif database.lower() == "ensembl":
            ref_type = 'ENSEMBL'
        elif database.lower() == "exontools":
            ref_type = 'EXONTOOLS'
        else:
            ref_type = 'CUSTOM'
        logging.info("Reference format type: {0:s}".format(ref_type))

        # PARSE REFERENCE FILE AND RETURN UNIQUE SEQUENCES
        TARGETS = makenewdir(name=os.path.join(tmpdir.path, "TARGETS"), fullname='TARGETS')
        parse_reference(refpath, database, organism, custom)
        annostats["targets"]["parsed"] = len(targref.bank)
        intargpath = os.path.join(TARGETS.path, 'inputref.fa')
        if not annotator.dryrun:
            with open(intargpath, 'w') as outfile:
                for gene in targref.bank:
                    outfile.write(">{0:s}\n{1:s}\n".format(gene, targref.bank[gene]['seq']))
        targclustpath = os.path.join(TARGETS.path, 'targclust.fa')
        outlog = os.path.join(log_dir.path, "target_clust.log")

        run_executor(executor(
            program="cd-hit-est",
            params=[intargpath, targclustpath, cluster, 50, mem, annotator.threads],
            conditions={"positive": [annotator.threads], "pathexists": [intargpath]},
            custom_arg_string=" > " + outlog + " 2>&1",
            quiet=True)
        )

        tclusters = []
        if not annotator.dryrun:
            with open(targclustpath, 'r') as outfile:
                for line in outfile:
                    if line.startswith(">"):
                        tclusters.append(line.strip()[1:])

        annostats["targets"]["clustered"] = len(tclusters)
        logging.info("Total {0:d} reference sequences remained after clusterization with {1:0.2f} identity threshold".format(len(tclusters), cluster))

        TARGETS.delete()
        TARGETS = makenewdir(name=os.path.join(tmpdir.path, "TARGETS"), fullname='TARGETS')

        outtargpath = os.path.join(TARGETS.path, 'targets.fa')
        targdb = os.path.join(TARGETS.path, 'targets')
        if not annotator.dryrun:
            with open(outtargpath, 'w') as outfile:
                for gene in sorted(tclusters, key=natural_sort):
                    outfile.write(">{0:s}\n{1:s}\n".format(gene, targref.bank[gene]['seq']))
            del tclusters

        # CREATE REFERENCE DATABASE
        outlog = os.path.join(log_dir.path, "target_makedb.log")
        run_executor(executor(
            program="makeblastdb",
            params=[outtargpath, targdb, target],
            conditions={"pathexists": [outtargpath]},
            custom_arg_string=" >> " + outlog + " 2>&1",
            quiet=True)
        )

        logging.info("Reference database is ready for annotation analysis")
        if annotator.debug:
            pdb.set_trace()

        logging.info("**********************************************************************")
        logging.info("STAGE 2. PREPARING QUERY SEQUENCES")
        annostats["queries"]["total"] = test_format(FileList, query)

        # FILTERING
        logging.debug("Perform filtration if references are provided")
        if filtering:
            FILTEMP = makenewdir(name=os.path.join(tmpdir.path, "FILTEMP"), fullname='FILTEMP')
            FileList = filtration(queries=FileList, filterpath=filtering, qtype=query, evalue=evalue, logpath=log_dir.path, tmppath=FILTEMP.path)

            annostats["queries"]["filtered"] = {}
            for f in FileList:
                count = 0
                libname = os.path.basename(f).split("_")[0].split(".")[0]
                ffasta = SeqIO(f)
                for seq in ffasta.read():
                    count += 1
                annostats["queries"]["filtered"][libname] = annostats["queries"]["total"][libname] - count
                del ffasta
            logging.debug("Finished filtration of query sequences: OK")
        if annotator.debug:
            pdb.set_trace()

        # CLUSTERING
        logging.info("Perform clusterization of query contigs with {0:0.2f} identity threshold".format(cluster))
        CLUSTERS = makenewdir(name=os.path.join(tmpdir.path, "CLUSTERS"), fullname='CLUSTERS')

        query_clusters = []
        annostats["queries"]["clustered"] = {}
        for f in FileList:
            outclustpath = os.path.join(CLUSTERS.path, os.path.basename(f))
            libname = os.path.basename(outclustpath).split("_")[0].split(".")[0]
            query_clusters.append(outclustpath)
            outlog = os.path.join(log_dir.path, "query_clust.log")
            run_executor(executor(
                program="cd-hit-est",
                params=[f, outclustpath, cluster, minlen, mem, annotator.threads],
                conditions={"positive": [annotator.threads], "pathexists": [f]},
                custom_arg_string=" > " + outlog + " 2>&1",
                quiet=True)
            )
            if not annotator.dryrun:
                os.remove(outclustpath + ".clstr")

            count = 0
            ffasta = SeqIO(f)
            for seq in ffasta.read():
                count += 1
            annostats["queries"]["clustered"][libname] = count
            logging.info("Total {0:d} query contigs remained in '{1:s}' library after clusterization".format(count, libname))
            del ffasta

        if filtering:
            FILTEMP.delete()
        logging.debug("Finished clusterization of query sequences: OK")

        if annotator.debug:
            pdb.set_trace()

        # ORF PREDICTION AND TRANSLATION
        if orf or translating:
            logging.info("Searching for long ORFs in each library")
            FileListBLAST = []
            if translating:
                logging.info("Each found ORF will be translated to peptide sequence")
            TASKS = []
            annostats["queries"]["ORFS"] = {}
            ORFS = makenewdir(name=os.path.join(tmpdir.path, "ORFS"), fullname='ORFS')

            # SEARCH ORFs within each file
            for fclust in query_clusters:
                TASKS.append(worker(search_orfs, [fclust, ORFS.path, translating, gencode, minlen, 0.5]))

            if TASKS:
                processes_requested = set_threads("ORF SEARCH", len(TASKS), annotator.threads)
                pool = create_pool(processes_requested)
                jobs = hard_worker(run_instance, TASKS, pool)
                close_pool(pool)

            if jobs:
                for result in jobs:
                    orfres.add_dict(result)
                del jobs
            else:
                logging.error("Multiprocessing error in ORF search")
                raise EXONtoolsError("Multiprocessing error in ORF search")

            if translating:
                query = "prot"
                annotat.translation()

            for lib in sorted(orfres.bank.keys(), key=natural_sort):
                FileListBLAST.append(orfres.bank[lib]["path"])
                logging.info("Total {0:d} ORFs were created in '{1:s}' library".format(orfres.bank[lib]["N_orf"], lib))
                annostats["queries"]["ORFS"][lib] = orfres.bank[lib]["N_orf"]

        else:
            FileListBLAST = query_clusters

        FileListBLAST = {os.path.basename(x).split("_")[0].split(".")[0]: x for x in FileListBLAST}

        logging.info("Query sequences are ready for annotation analysis")

        if annotator.debug:
            pdb.set_trace()

        logging.info("**********************************************************************")
        logging.info("STAGE 3. BLAST ALIGNMENT")

        outfmtdict = {}

        BLAST = makenewdir(name=os.path.join(tmpdir.path, "BLAST"), fullname='BLAST')
        outlog = os.path.join(log_dir.path, "blast.log")
        if query == "prot" and target == "prot":
            blastprog = "blastp"
        elif query == "nucl" and target == "nucl":
            blastprog = "blastn"
        elif query == "nucl" and target == "prot":
            blastprog = "blastx"
        elif query == "prot" and target == "nucl":
            blastprog = "tblastn"
        else:
            logging.error("Unknown BLAST algorithm")
            raise EXONtoolsError("Unknown BLAST algorithm")

        for lib in sorted(FileListBLAST.keys(), key=natural_sort):
            outfmtdict[lib] = blastalign(lib, FileListBLAST[lib], targdb, blastprog, evalue, gencode, outlog, BLAST.path)

        logging.debug("BLAST alignment completed: OK")

        if annotator.debug:
            pdb.set_trace()

        logging.info("**********************************************************************")
        logging.info("STAGE 4. RUNNING ANNOTATION ANALYSIS")

        annotat.setblast(blastprog)

        TASKS = []
        for lib in sorted(outfmtdict.keys(), key=natural_sort):
            try:
                if chimeras:
                    chlib = chimeras[lib]
                else:
                    chlib = False
                if scaffolds:
                    scaflib = scaffolds[lib]
                else:
                    scaflib = False
                if isoforms:
                    isolib = isoforms[lib]
                else:
                    isolib = False
                TASKS.append(worker(parse_blast, [lib, outfmtdict[lib], similarity, overlap, scaflib, chlib, isolib]))
            except KeyError:
                logging.error("All names of optional files must be the same as the names of annotation libraries")
                raise EXONtoolsError("Optional file name and annotation library name do not match!")

        if TASKS:
            processes_requested = set_threads("ANNOTATION", len(TASKS), annotator.threads)
            pool = create_pool(processes_requested)
            jobs = hard_worker(run_instance, TASKS, pool)
            close_pool(pool)

        if jobs:
            for result in jobs:
                annotres.add_dict(result)
            del jobs

            annotres.validate()

        else:
            logging.error("Multiprocessing error in BLAST analysis")
            raise EXONtoolsError("Multiprocessing error in BLAST analysis")

        # annotout = makenewdir(name="ANNOTATION", fullname='ANNOTATION')
        # tandemout = makenewdir(name="REPEATS", fullname='REPEATS')
        # unknownout = makenewdir(name="UNKNOWN", fullname='UNKNOWN')
        # longorfout = makenewdir(name=os.path.join(unknownout.path,"ORFS"), fullname='long ORFS')

        logging.debug("Annotation results successfully parsed: OK")

        if annotator.debug:
            pdb.set_trace()

        logging.info("**********************************************************************")
        logging.info("STAGE 5. SAVING ANNOTATION RESULTS")
        TASKS = []
        annostats["queries"]["BLASTED"] = {}
        annostats["queries"]["ANNOTATED"] = {}
        annostats["queries"]["UNIQUE"] = {}
        annostats["queries"]["UNKNOWN"] = {}
        annostats["queries"]["UNKNORFS"] = {}
        annostats["queries"]["REPEATS"] = {}
        annostats["queries"]["SCAFFOLDS"] = {}
        annostats["queries"]["ISOFORMS"] = {}
        annostats["queries"]["CHIMERAS"] = {}
        libstats = {}

        for fqc in sorted(query_clusters, key=natural_sort):
            libname = os.path.basename(fqc).split("_")[0].split(".")[0]
            TASKS.append(worker(finalize_annotations,[libname,fqc,norename,translating, oneway,gencode, minlen, 0.5]))
            annostats["queries"]["BLASTED"][libname] = len(annotres.bank[libname])
            logging.info("Total {0:d} queries were annotated in '{1:s}' library".format(len(annotres.bank[libname]), libname))

        if TASKS:
            processes_requested = set_threads("FINALIZE", len(TASKS), annotator.threads)
            pool = create_pool(processes_requested)
            jobs = hard_worker(run_instance, TASKS, pool)
            close_pool(pool)

        if jobs:
            for result in jobs:
                libname = list(result.keys())[0]
                annostats["queries"]["ANNOTATED"][libname] = result[libname][0]
                annostats["queries"]["SCAFFOLDS"][libname] = result[libname][1]
                annostats["queries"]["ISOFORMS"][libname] = result[libname][2]
                annostats["queries"]["CHIMERAS"][libname] = result[libname][3]
                annostats["queries"]["UNKNOWN"][libname] = result[libname][4]
                annostats["queries"]["UNKNORFS"][libname] = result[libname][5]
                annostats["queries"]["REPEATS"][libname] = result[libname][6]
                annostats["queries"]["UNIQUE"][libname] = result[libname][7]
                libstats[libname]=[]
            del jobs
        else:
            logging.error("Multiprocessing error in parsing ANNOTATION results")
            raise EXONtoolsError("Multiprocessing error in parsing ANNOTATION results")

        if annotator.debug:
            pdb.set_trace()

        if annotator.stats:
            logging.info("**********************************************************************")
            logging.info("STAGE 6. COLLECTING ANNOTATION STATS")
            for lib in libstats:
                libstats[lib].append(annostats["queries"]["total"][lib])
                if "filtered" in annostats["queries"]:
                    libstats[lib].append(annostats["queries"]["filtered"][lib])
                libstats[lib].append(annostats["queries"]["clustered"][lib])
                if "ORFS" in annostats["queries"]:
                    libstats[lib].append(annostats["queries"]["ORFS"][lib])
                libstats[lib].append(annostats["queries"]["BLASTED"][lib])
                libstats[lib].append(annostats["queries"]["ANNOTATED"][lib])
                libstats[lib].append(annostats["queries"]["UNIQUE"][lib])
                libstats[lib].append(annostats["queries"]["SCAFFOLDS"][lib])
                libstats[lib].append(annostats["queries"]["ISOFORMS"][lib])
                libstats[lib].append(annostats["queries"]["CHIMERAS"][lib])
                libstats[lib].append(annostats["queries"]["UNKNOWN"][lib])
                libstats[lib].append(annostats["queries"]["UNKNORFS"][lib])
                libstats[lib].append(annostats["queries"]["REPEATS"][lib])
            statpath = os.path.join(output.path, 'annotation_stats.csv')
            logging.info("All annotation stats will be saved to 'annotation_stats.csv' file")
            header = ["No","RRL","TOTAL","FILTERED","CLUSTERED","ORFs","ALIGNED","ANNOTATED","UNIQUE","SCAFFOLDS","ISOFORMS","CHIMERAS","UNKNOWN","UNKNOWN_ORFs","UNKNOWN_STRs"]
            if "ORFS" not in annostats["queries"]:
                header.remove("ORFs")
            if "filtered" not in annostats["queries"]:
                header.remove("FILTERED")
            if not annotator.dryrun:
                with open(statpath, 'w') as statfile:
                    csv_writer = csv.writer(statfile)
                    csv_writer.writerow(header)
                    for i,lib in enumerate(sorted(libstats.keys(), key=natural_sort)):
                        csv_writer.writerow([i+1,lib]+libstats[lib])
            logging.debug("Annotation stats were successfully written to the file: OK")

        if not annotator.keeptmp:
            tmpdir.delete()


def finalize_annotations(lib,inpath,norename,translating,oneway,genecode,minlen,difforf):
    """Writes annotation results to files"""

    origfasta = SeqIO(inpath)
    annotoutpath = os.path.join(output.path,lib+"_annotation.gff")
    annotcontigpath = os.path.join(output.path,lib+"_annotation.fasta")
    unknowncontigpath = os.path.join(output.path,lib+"_unknown.fasta")
    STRpath = os.path.join(output.path,lib+"_STR.fasta")
    longorfpath = os.path.join(output.path,lib+"_unknownorfs.fasta")

    annotID = 0
    genes = []
    scafcount = 0
    isocount = 0
    chimercount = 0
    unknowncount = 0
    uknownorf = 0
    repcount = 0
    annotated = {}

    if not annotator.dryrun:
        with open(annotoutpath,'a') as outfile:
            outfile.write("##gff-version 3\n")
            outfile.write("#!gff-spec-version 1.23\n")
            outfile.write("#!EXONtools pipeline\n")
            outfile.write("#!EXONtools command annotate_contigs\n")
            outfile.write("#!alignment processor NCBI {0:s}\n".format(annotat.blastprog))
            outfile.write("#!input assembly {0:s}\n".format(os.path.basename(inpath)))
            outfile.write("#!annotation reference {0:s}\n".format(os.path.basename(targref.path)))
            outfile.write("#!IMPORTANT: gene coordinates refer to annotated regions in each query contig\n".format(os.path.basename(targref.path)))

        for contig in origfasta.read():
            results = []
            contigname = contig.name.split()[0]
            if orfres.bank:
                for orfseq in orfres.bank[lib]['contigs'][contigname]:
                    try:
                        resannot = annotres.bank[lib][contigname+"_"+orfseq["name"]]
                        resannot.addorf(orfseq)
                        results.append(resannot)
                    except KeyError:
                        pass
            else:
                try:
                    results.append(annotres.bank[lib][contigname])
                except KeyError:
                        pass
            if results:
                results.sort(key=lambda x: x.score)
                annotated[contigname] = results[0]
                with open(annotoutpath,'a') as outfile:
                    outfile.write("##sequence-region\t{0:s}\t1\t{1:d}\n".format(contigname,len(contig.seq)))
            else:
                orflist = findorfs(seq=contig.seq, code=genecode, minlen=minlen, diff=difforf)
                if orflist:
                    uknownorf +=1
                    orflist.sort(key = lambda x: x.length, reverse=True)
                    orfseq = orflist[0]
                    with open(longorfpath,'a') as outfile:
                        outfile.write(">{0:s}_ORF_{1:d}_{2:d}_[{3:d}]\n{4:s}\n".format(contigname,orfseq.start,orfseq.end,orfseq.frame,orfseq.seq))

                repeats = microsat(contig.seq, minlen=2, maxlen=21, minrep=3, minseqlen=9, mismatch=0.1)
                if repeats:
                    repcount+=1
                    with open(STRpath,'a') as outfile:
                        for r in repeats:
                            outfile.write("{0:s}\t{1:d}\t{2:d}\t{3:d}\t{4:s}\t{5:s}\n".format(contigname,r[0],r[1],r[2],r[3],r[4]))

                with open(unknowncontigpath,'a') as outfile:
                    outfile.write("{0:s}\n{1:s}\n".format(contig.name,contig.seq))

        if annotat.is_translated:
            trans_adjust = 3
        else:
            trans_adjust = 1

        for contig in origfasta.read():
            contigname = contig.name.split()[0]
            contlength = len(contig.seq)
            try:
                result = annotated[contigname]
                if result.orf:
                    if result.orf["frame"] > 0:
                        orfframe = "+"
                        cstart = result.orf["start"] + (result.start - 1)*trans_adjust
                        cend = result.orf["start"] - 1 + result.end *trans_adjust
                        cframe = result.strand
                    else:
                        orfframe = "-"
                        cend = result.orf["end"] - (result.start-1) *trans_adjust
                        cstart = result.orf["end"]  + 1 - result.end *trans_adjust
                        cframe = plustominus(result.strand)

                else:
                    cstart = result.start
                    cend = result.end
                    cframe = result.strand

                annotID +=1
                genes.append(targref.bank[result.target]["gene"])
                with open(annotcontigpath,'a') as outfile:
                    if not norename:
                        extrapart = ""
                        if result.scaffolds:
                            extrapart = extrapart + ";scaffold segment"
                        if result.isoforms:
                            extrapart = extrapart + ";isoform"
                        if result.chimeric:
                            extrapart = extrapart + ";chimeric"

                        if oneway and cframe == "-":
                            outfile.write(">{0:s}\ttarget={1:s};name={2:s};gene={3:s};organism={4:s};description={5:s}{6:s}\n{7:s}\n".format(contig.name,targref.bank[result.target]["id"],targref.bank[result.target]["gene"],targref.bank[result.target]["name"], targref.bank[result.target]["organism"],targref.bank[result.target]["description"],extrapart, DNArevcomp(contig.seq)))
                        else:
                            outfile.write(">{0:s}\ttarget={1:s};name={2:s};gene={3:s};organism={4:s};description={5:s}{6:s}\n{7:s}\n".format(contig.name,targref.bank[result.target]["id"],targref.bank[result.target]["gene"],targref.bank[result.target]["name"], targref.bank[result.target]["organism"],targref.bank[result.target]["description"],extrapart,contig.seq))
                    else:
                        if oneway and cframe == "-":
                            outfile.write(">{0:s}\n{1:s}\n".format(contig.name,DNArevcomp(contig.seq)))
                        else:
                            outfile.write(">{0:s}\n{1:s}\n".format(contig.name,contig.seq))
                with open(annotoutpath,'a') as outfile:
                    annotstring = "ID=gene{0:d};Target={1:s} {2:d} {3:d};Name={4:s};gene={7:s};organism={5:s};description={6:s}".format(annotID,targref.bank[result.target]["id"],result.tstart,result.tend,targref.bank[result.target]["gene"],targref.bank[result.target]["organism"],targref.bank[result.target]["description"],targref.bank[result.target]["name"])

                    cleanedscafs = []
                    if result.scaffolds:
                        cleanedscafs = [re.sub("_ORF\d+","",x) for x in result.scaffolds]
                        cleanedscafs = list(set([x for x in cleanedscafs if x != contigname]))
                        if cleanedscafs:
                            annotstring = annotstring + ";scaffolds="+",".join(cleanedscafs)
                            scafcount +=1
                    if result.isoforms:
                        cleanediso = [re.sub("_ORF\d+","",x) for x in result.isoforms]
                        cleanediso = list(set([x for x in cleanediso if x != contigname and x not in cleanedscafs]))
                        if cleanediso:
                            annotstring = annotstring + ";isoforms="+",".join(cleanediso)
                            isocount+=1
                    if result.chimeric:
                        annotstring = annotstring + ";is_chimeric=true"
                        chimercount+=1

                    if oneway and cframe == "-":
                        outfile.write("{0:s}\t{1:s}\tgene\t{2:d}\t{3:d}\t{4:.2e}\t{5:s}\t.\t{6:s}\n".format(contigname,result.source,contlength-cend+1,contlength-cstart+1,Decimal(result.score),"+",annotstring))
                        if result.orf:
                            outfile.write("{0:s}\t{1:s}\tCDS\t{2:d}\t{3:d}\t.\t{4:s}\t{5:d}\tID=CDS{6:d};Parent=gene{6:d}\n".format(contigname,result.source,contlength-result.orf["end"]+1,contlength-result.orf["start"]+1,plustominus(orfframe),abs(result.orf['frame'])-1,annotID))
                    else:
                        outfile.write("{0:s}\t{1:s}\tgene\t{2:d}\t{3:d}\t{4:.2e}\t{5:s}\t.\t{6:s}\n".format(contigname,result.source,cstart,cend,Decimal(result.score),cframe,annotstring))
                        if result.orf:
                            outfile.write("{0:s}\t{1:s}\tCDS\t{2:d}\t{3:d}\t.\t{4:s}\t{5:d}\tID=CDS{6:d};Parent=gene{6:d}\n".format(contigname,result.source,result.orf["start"],result.orf["end"],orfframe,abs(result.orf['frame'])-1,annotID))
            except KeyError:
                unknowncount +=1
    genes = [x.split("_")[0] for x in genes]
    genes = list(set(genes))

    return {lib:(annotID,scafcount,isocount,chimercount,unknowncount,uknownorf,repcount,len(genes))}


def parse_blast(libname, blastfile, similarity, overlap, scafpath=False, chimpath=False, isopath=False):
    """Performs annotation of a particular library based on blast results"""

    logging.info("Parsing annotation results for '{0:s}' library".format(libname))

    # READ BLAST OUTPUT
    blastout = BLAST6(blastfile)

    if annotat.blastprog == 'blastn':
        scafover=24
    else:
        scafover=8

    targets = {}
    queries = {}
    ANNOTATION = {}

    # PARSING EXTRA FILES
    scaffolds = parse_extra(scafpath)
    isoforms = parse_extra(isopath)
    chimeras = parse_extra(chimpath)

    if chimeras:
        # make a list with all chimeric contigs
        chimeras = list(set(list(chimeras.keys()) + sum(chimeras.values(), [])))

    for hitline in blastout.read():
        try:
            targets[hitline.target].append(hitline)
        except KeyError:
            targets[hitline.target] = [hitline]
        try:
            queries[hitline.query].append(hitline)
        except KeyError:
            queries[hitline.query] = [hitline]

    for contig in queries:
        # try to annotate each contig in queries according to all possible scenarios

        # PARSE QUERIES
        if len(queries[contig]) > 1:
            # filter queries

            filtered_queries = []
            for ref in queries[contig]:
                if ref.length >= overlap and ref.identity >= similarity and not [x for x in filtered_queries if x.target == ref.target]:
                    filtered_queries.append(ref)
                else:
                    # remove filtered matches from target dictionary
                    try:
                        targets[ref.target] = [x for x in targets[ref.target] if x.id != ref.id]
                        if not targets[ref.target]:
                            del targets[ref.target]
                    except KeyError:
                        pass

            # Select the best annotation among filtered querires
            # if the number of filtered contigs is greater 1 - test chimeras
            if filtered_queries:
                selected_query = filtered_queries[0]
                # check misassemblies
                if test_chimeras(filtered_queries, chimeras):
                    selected_query.missasembled()

                queries[contig] = selected_query
                for ref in filtered_queries[1:]:
                    # remove remaining matches from target dictionary
                    try:
                        targets[ref.target] = [x for x in targets[ref.target] if x.id != ref.id]
                        if not targets[ref.target]:
                            del targets[ref.target]
                    except KeyError:
                        pass
            else:
                continue
        elif len(queries[contig]) == 1:
            queries[contig] = queries[contig][0]
        else:
            continue

        # PARSE TARGETS
        qtarget = queries[contig].target
        try:
            # Is it just a single annotation?
            if len(targets[qtarget]) == 1:
                if queries[contig].length >= overlap and queries[contig].identity >= similarity:
                    ANNOTATION[contig] = annotat(queries[contig])
                targets[qtarget] = [x for x in targets[qtarget] if x.id != queries[contig].id]
                if not targets[qtarget]:
                    del targets[qtarget]

            # if a reference target has multiple query hits
            else:

                # filter refs
                filtered_refs = []
                for ref in targets[qtarget]:
                    if ref.length >= overlap and ref.identity >= similarity and not [x for x in filtered_refs if x.query == ref.query]:
                        filtered_refs.append(ref)

                # CHECK ALL POSSIBLE SCENARIOS:
                # delete target if all refs were filtered
                if not filtered_refs:
                    del targets[qtarget]

                # if the remaining ref is our query contig
                elif len(filtered_refs) == 1 and queries[contig].id == filtered_refs[0].id:
                    ANNOTATION[contig] = annotat(queries[contig])
                    del targets[qtarget]

                # if selected target has multiple refs including query contig
                elif [x for x in filtered_refs if x.query == contig]:
                    ANNOTATION[contig] = annotat(queries[contig])

                    # check scaffolding
                    filtered_refs = [x for x in filtered_refs if x.query != contig]
                    scaflist = test_scaffolds(filtered_refs, queries[contig], scaffolds, maxoverlap = scafover)
                    if scaflist:
                        [x.add_scaf(contig) for x in filtered_refs if x.query in scaflist]
                        [ANNOTATION[contig].add_scaf(x) for x in scaflist]
                        isolist = [x.query for x in filtered_refs if x.query not in scaflist]
                    else:
                        isolist = [x.query for x in filtered_refs]

                    # check isoforms
                    if isoforms:
                        try:
                            isochecks = isoforms[contig]
                            isolist = [x for x in isolist if x in isochecks]
                        except KeyError:
                            pass
                    [x.add_iso(contig) for x in filtered_refs if x.query in isolist]
                    [ANNOTATION[contig].add_iso(x) for x in isolist]
                    targets[qtarget] = filtered_refs

                # update target refs if they do not contain our query contig
                else:
                    targets[qtarget] = filtered_refs

        # if the target instance has been already parsed and deleted
        except KeyError:
            continue

    if targets:
        logging.error("Not all query annotations were analyzed. Annotator program bug.")
        logging.error(list(targets.keys()))
        raise EXONtoolsError("Annotator program bug.")

    return {libname: ANNOTATION}


class annotat(object):
    """A single annotation instance"""

    source = "EXONtools"
    is_translated = False
    blastprog = None
    dryrun = False
    reference = None

    def __init__(self, blastline):
        self.sequence = blastline.query
        self.feature = None
        self.start = min(blastline.qstart, blastline.qend)
        self.end = max(blastline.qstart, blastline.qend)
        self.tstart = min(blastline.tstart, blastline.tend)
        self.tend = max(blastline.tstart, blastline.tend)
        self.score = blastline.evalue
        self.orf = None
        if blastline.qstart < blastline.qend and blastline.tstart < blastline.tend:
            self.strand = "+"
        elif blastline.qstart > blastline.qend or blastline.tstart > blastline.tend:
            self.strand = "-"
        else:
            self.strand = "."
        self.target = blastline.target
        self.scaffolds = blastline.scaffolds
        self.isoforms = blastline.isoforms
        self.chimeric = blastline.chimeric

    @classmethod
    def setblast(cls, prog):
        if isinstance(prog, str):
            logging.debug("'{0:s}' BLAST program has been assigned to 'annotation' class".format(prog))
            cls.blastprog = prog
        else:
            logging.error("Input type error in annotat.setblast method")
            raise EXONtoolsError("Input type error in annotat.setblast method")

    @classmethod
    def setsource(cls, name):
        if isinstance(name, str):
            logging.debug("'{0:s}' program has been assigned to source tag".format(name))
            cls.source = name
        else:
            logging.error("Input type error in annotat.setsource method")
            raise EXONtoolsError("Input type error in annotat.setsource method")

    @classmethod
    def translation(cls):
        cls.is_translated = True

    @classmethod
    def setdry(cls):
        cls.dryrun = True

    def addorf(self,orf):
        self.orf = orf

    def missasembled(self):
        self.chimeric = True

    def add_iso(self, name):
        self.isoforms.append(name)

    def add_scaf(self, name):
        self.scaffolds.append(name)


def test_scaffolds(blastlist, query, scafdict, maxoverlap=25):
    """test if there are any scaffolds"""

    scaffolds = []
    qset = set(range(query.tstart, query.tend + 1))

    for ref in blastlist:
        refset = set(range(ref.tstart, ref.tend + 1))
        if len(qset & refset) <= maxoverlap:
            scaffolds.append(ref.query)

    if scafdict:
        try:
            scafs = scafdict[query.query]
            scaffolds = [x for x in scaffolds if x in scafs]
        except KeyError:
            pass

    return scaffolds


def test_chimeras(blastlist, chimeras=False, maxoverlap=10):
    """test if query represents chimeric sequence"""

    if len(blastlist) < 2:
        return False
    elif chimeras and blastlist[0].query not in chimeras:
        return False
    else:
        for i, ref in enumerate(blastlist):
            refset = set(range(min(ref.qstart, ref.qend), max(ref.qstart, ref.qend) + 1))
            for check in blastlist[i + 1:]:
                checkset = set(range(min(check.qstart, check.qend), max(check.qstart, check.qend) + 1))
                if len(refset & checkset) <= maxoverlap:
                    return True


def parse_extra(inpath):
    """Parsing scaffolds, isoforms and chimeric files"""

    if inpath and os.path.exists(inpath):
        outdict = {}
        with open(inpath, 'r') as infile:
            for line in infile:
                x = line.strip().split("\t")
                contig = x[0].strip()
                outdict[contig] = sorted([y.strip() for y in x[1].split(",")], key=natural_sort)
        return outdict
    else:
        return False


def blastalign(libname, infile, refdb, blastprog, evalue, gencode, logpath, tmppath):
    """Align queries against references with BLAST"""

    blastout = os.path.join(tmppath, libname + ".outfmt6")
    logging.info("Running '{0:s}' analysis for '{1:s}' library".format(blastprog, libname))
    if blastprog == 'tblastn' or blastprog == "blastx":
        blastpar = "-query_gencode " + str(gencode)
    else:
        blastpar = ""

    run_executor(executor(
        program=blastprog,
        params=[refdb, infile, blastout, 6, str(evalue), annotator.threads, blastpar],
        conditions={"positive": [annotator.threads, gencode], "pathexists": [infile]},
        custom_arg_string=annotator.extra + " >> " + logpath + " 2>&1",
        quiet=True)
    )
    return blastout


def search_orfs(inpath, outpath, translation, genecode, minlen, difforf):
    """Provides the dictionary with all found ORFs"""

    if translation:
        def orfrecord(outfile, contig, count, orfseq, code):
            outfile.write(">{0:s}_ORF{1:d}\n{2:s}\n".format(contig, count, translate(orfseq.seq, 1, code)))
    else:
        def orfrecord(outfile, contig, count, orfseq, code):
            outfile.write(">{0:s}_ORF{1:d}\n{2:s}\n".format(contig, count, orfseq.seq))

    if os.path.exists(inpath) and os.path.exists(outpath):
        libname = os.path.basename(inpath).split("_")[0].split(".")[0]
        logging.debug("Searching long ORFs in '{0:s}' library".format(libname))
        infile = SeqIO(inpath)
        filepath = os.path.join(outpath, libname + "_orfs.fasta")

        ORFdict = {libname: {"path": filepath, "N_orf": 0, "contigs": {}}}

        with open(filepath, 'w') as orffile:
            for inseq in infile.read():
                orflist = findorfs(seq=inseq.seq, code=genecode, minlen=minlen, diff=difforf)
                contig = inseq.name.split()[0]
                ORFdict[libname]["contigs"][contig]=[]
                count = 0
                if orflist:
                    orflist.sort(reverse=True, key=lambda x: x.length)
                    for orfseq in orflist:
                        count += 1
                        orfrecord(orffile, contig, count, orfseq, genecode)
                        ORFdict[libname]["contigs"][contig].append({"name":"ORF"+str(count),"length":orfseq.length,"start":orfseq.start,"end":orfseq.end,"frame":orfseq.frame})
                ORFdict[libname]["N_orf"] += count
        return ORFdict
    elif annotator.dryrun:
        libname = os.path.basename(inpath).split("_")[0].split(".")[0]
        filepath = os.path.join(outpath, libname + "_orfs.fasta")
        ORFdict = {libname: {"path": filepath, "N_orf": 0, "contigs": {}}}
        return ORFdict
    else:
        logging.error("IO path for ORF search does not exist")
        raise EXONtoolsError("Search ORF IO path error")


def filtration(queries, filterpath, qtype, evalue, logpath, tmppath):
    """Filter query contigs based on provided references"""

    getinput.format(['.fasta', '.pep'])

    if os.path.exists(filterpath):
        logging.info("Filtering query contigs based on provided references")
        getinput.setdry()
        filtlist = getinput(filterpath).files
        getinput.undry()
    else:
        logging.error("The provided path for filtering references does not exist")
        raise EXONtoolsError("Filtering input path error")

    if not filtlist:
        logging.error("The provided path for filtering references does not contain supported files")
        raise EXONtoolsError("Filtering input path error")

    FILTERED = makenewdir(name="FILTERED", fullname='FILTERED')

    outlog = os.path.join(logpath, "query_filtration.log")

    logging.debug("Create blast databases and choose the blast algorithm")

    if annotator.debug:
        pdb.set_trace()

    filterdbs = []
    for filtpath in filtlist:
        filtrefname = os.path.splitext(os.path.basename(filtpath))[0].split("_")[0]
        dbnamepath = os.path.join(tmppath, filtrefname)

        if filtpath.endswith(".fasta"):
            dbtype = 'nucl'
        elif filtpath.endswith(".pep"):
            dbtype = 'prot'
        else:
            logging.warning("Unknown reference file format. Skipping filtering against '{0:s}'.".format(os.path.basename(filtpath)))

        run_executor(executor(
            program="makeblastdb",
            params=[filtpath, dbnamepath, dbtype],
            conditions={"pathexists": [filtpath]},
            custom_arg_string=" >> " + outlog + " 2>&1",
            quiet=True)
        )

        logging.debug("Choose alignment algorithm for filtration")
        if qtype == 'nucl' and dbtype == 'nucl':
            blprog = 'blastn'
        elif qtype == 'nucl' and dbtype == 'prot':
            blprog = 'blastx'
        elif qtype == 'prot' and dbtype == 'nucl':
            blprog == 'tblatn'
        elif qtype == 'prot' and dbtype == 'prot':
            blprog == 'blastp'
        else:
            raise EXONtoolsError("Uknown blast algorithm for filtering")

        filterdbs.append((dbnamepath, blprog))

    logging.debug("Perform filtration alignments")
    if annotator.debug:
        pdb.set_trace()

    newfiles = []
    for f in queries:

        libname = os.path.splitext(os.path.basename(f))[0].split("_")[0]
        filtpath = os.path.join(tmppath, os.path.basename(f))
        newfiles.append(filtpath)

        infasta = SeqIO(f)
        inseqs = []
        for inseq in infasta.read():
            inseqs.append(inseq)
        del infasta

        totalfilter = []
        for fdb in filterdbs:
            filtresult = fdb[0] + '.outfmt6'
            filteredpath = os.path.join(FILTERED.path, libname + "_" + os.path.basename(fdb[0]) + "_filtered.fa")
            run_executor(executor(
                program=fdb[1],
                params=[fdb[0], f, filtresult, 6, str(evalue), annotator.threads, ""],
                conditions={"positive": [annotator.threads], "pathexists": [f]},
                custom_arg_string=annotator.extra + " >> " + outlog + " 2>&1",
                quiet=True)
            )

            blastres = BLAST6(filtresult)
            namestofilter = [x.query for x in blastres.read()]
            namestofilter = list(set(namestofilter))

            # blastres.delete()
            del blastres

            # save filtered contigs
            filterlist = []
            if not annotator.dryrun:
                with open(filteredpath, 'w') as outfile:
                    for i, inseq in enumerate(inseqs):
                        if inseq.name.split()[0] in namestofilter:
                            filterlist.append(i)
                            outfile.write(">{0:s}\n{1:s}\n".format(inseq.name, inseq.seq))

            totalfilter = totalfilter + filterlist

        logging.debug("Parsing filtering results for '{0:s}' library".format(libname))

        # delete filtered contigs
        totalfilter = list(set(totalfilter))
        totalfilter.sort(reverse=True)
        logging.info("{0:d} query contigs were filtered in '{1:s}' library".format(len(totalfilter), libname))
        for i in totalfilter:
            del inseqs[i]

        # save remained contigs
        if not annotator.dryrun:
            with open(filtpath, 'w') as outfile:
                for inseq in sorted(inseqs, key=lambda x: natural_sort(x.name)):
                    outfile.write(">{0:s}\n{1:s}\n".format(inseq.name, inseq.seq))
            del inseqs

    return newfiles


def test_format(inpath, stype):
    """Function that tests nucleotide input files"""

    results = {}
    cutnucls = set("ACGT")
    nucls = set("ACGTNRYSWKMBDHV")
    prots = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*")

    if stype == "nucl":
        def test(seq, path):
            if set(seq) - nucls:
                logging.error("'{0:s}'".format(seq))
                logging.error("Selected sequence type 'nucl' is not confirmed in '{0:s}'".format(os.path.basename(path)))
                logging.error("Please choose the correct sequence type of query and target sequences")
                raise EXONtoolsError("Sequence type error")
    elif stype == "prot":
        def test(seq, path):
            if len(seq)>100 and (set(seq) - prots or not set(seq) - cutnucls):
                logging.error("'{0:s}'".format(seq))
                logging.error("Selected sequence type 'prot' is not confirmed in '{0:s}'".format(os.path.basename(path)))
                logging.error("Please choose the correct sequence type of query and target sequences")
                raise EXONtoolsError("Sequence type error")
    else:
        logging.error("Unknown format provided for input file test")
        raise EXONtoolsError("Unknown file format")

    for x in inpath:
        testin = SeqIO(x)
        count = 0
        for seq in testin.read():
            count += 1
            test(seq.seq, x)
        libid = os.path.basename(x).split("_")[0].split(".")[0]
        results[libid] = count
        logging.info("Library ID: {0:s} | Seq type: {1:s} | Num seqs: {2:d}".format(libid, stype, count))
    return results


class targref(object):
    """Class to store all targets"""

    path = None
    bank = {}

    @classmethod
    def add_dict(cls,tdict):
        cls.bank = tdict

    @classmethod
    def add_path(cls,inpath):
        cls.path = inpath


class annotres(object):
    """Class to store all annotation results"""

    bank = {}

    @classmethod
    def add_dict(cls,tdict):
        cls.bank.update(tdict)

    @classmethod
    def validate(cls):
        logging.info("Running the final validation for scaffolds and isoforms")
        #test scaffolding contigs:
        for lib in cls.bank:
            check = set(list(cls.bank[lib].keys()))
            for contig in cls.bank[lib]:
                if set(cls.bank[lib][contig].scaffolds)-check:
                    logging.warning("Scaffold segment '{0:s}' of '{1:s}' contig is not present in the annotation of '{2:s}' library".format(' '.join(list(set(cls.bank[lib][contig].scaffolds)-check)),cls.bank[lib][contig].sequence,lib))
                    # raise EXONtoolsError("Annotation validation error")
                if set(cls.bank[lib][contig].isoforms)-check:
                    logging.warning("Isoform '{0:s}' of '{1:s}' contig is not present in the annotation of '{2:s}' library".format(' '.join(list(set(cls.bank[lib][contig].isoforms)-check)), cls.bank[lib][contig].sequence,lib))
                    # raise EXONtoolsError("Annotation validation error")


class orfres(object):
    """Class to store all ORF results"""

    bank = {}

    @classmethod
    def add_dict(cls,tdict):
        cls.bank.update(tdict)


def parse_reference(reference, database, organism, custom):
    """Parse reference based on predefined or user-provided grep patterns.
    The algorithm takes the longest sequence (scaffold) for each protein
    Returns the dictionary with proteins and captured sequence parameters
    """

    logging.debug("Parsing the reference database")

    # DEFINE WHICH GREP PATTERN TO USE
    if database.lower() == "ensembl":
        ref_code_pattern = "^(ENS\S+)"
        gene_id_pattern = "gene:(ENS\S+)"
        ref_name_pattern = "gene_symbol:(\S+)"
        ref_organism_pattern = None
        description_pattern = "description:([^\[]+) \["

    elif database.lower() == "swissprot":
        ref_code_pattern = "^sp\|(\w+)\|\w+"
        gene_id_pattern = "^sp\|\w+\|(\w+)"
        ref_name_pattern = "GN=(\w+)"
        ref_organism_pattern = "OS=([^=]+)\s+\w+="
        description_pattern = "^sp\|\w+\|\w+\s+([^=]+)\s+OS="

    elif database.lower() == "exontools":
        ref_code_pattern = "^([^\s]+).*$"
        gene_id_pattern = "^[^\s]+\s+.+name=([^;]+);.*$"
        ref_name_pattern = "^[^\s]+\s+.+gene=([^;]+);.*$"
        ref_organism_pattern = "^[^\s]+\s+.+organism=([^;]+);.*$"
        description_pattern = "^[^\s]+\s+.+description=([^;]+);.*$"

    elif database.lower() == "custom":
        ref_code_pattern = "".join(re.split("\)|\(", custom[::-1], 2))[::-1]
        gene_id_pattern = "".join(re.split("\)|\(", custom, 2))
        ref_name_pattern = None
        ref_organism_pattern = None
        description_pattern = None

    else:
        raise EXONtoolsError("Unknown reference database parsing format")

    # PARSING THE FILE WITH REFERENCE

    target_dict = {}                    # dictionary for storing selected targets
    targcount = 0
    reffile = SeqIO(reference)

    for refseq in reffile.read():

        if re.search("pseudogene", refseq.name.lower()):
            continue

        try:
            seq_id = re.search(ref_code_pattern, refseq.name).group(1).strip()             # get target sequence ID
            gene_id = re.search(gene_id_pattern, refseq.name).group(1).strip()              # get GENE ID
            seq_id = re.sub(",|;|=","",seq_id)
            gene_id = re.sub(",|;|=","",gene_id)

            if ref_name_pattern:
                try:
                    ref_name = re.search(ref_name_pattern, refseq.name).group(1)   # get reference protein ID
                    ref_name = re.sub(",|;|=","",ref_name)
                except AttributeError:
                    ref_name = 'NA'
            else:
                ref_name = 'NA'
            if organism:
                org_name = organism
            elif ref_organism_pattern:
                try:
                    org_name = re.search(ref_organism_pattern, refseq.name).group(1)   # get reference organism ID
                except AttributeError:
                    org_name = 'NA'
            else:
                org_name = 'NA'
            if description_pattern:
                try:
                    seq_info = re.search(description_pattern, refseq.name).group(1)  # get additional info
                    seq_info = re.sub(",|;|=","",seq_info)
                except AttributeError:
                    seq_info = 'NA'
            else:
                seq_info = 'NA'

        # THIS JUST CHECKS FOR SOME POSSIBLE BUGS - PROVIDED GREP PATTERNS MUST MATCH THE DATA
        except AttributeError as e:
            logging.error(refseq.name)
            logging.error(e)
            logging.error("Selected grep pattern can't match the reference gene.")
            logging.error(" Please verify that you define the correct type of reference database in your input options.")
            raise EXONtoolsError("EXONtools CRITICAL ERROR")

        targcount +=1
        target_dict["target" + str(targcount)] = {"id": seq_id, "seq": refseq.seq, "length": len(refseq.seq), "name": ref_name, "organism": org_name, "description": seq_info, "gene": gene_id}


    #     try:
    #         # THIS JUST CHECKS FOR SOME POSSIBLE BUGS - EACH SELECTED GENE TARGET MUST BE UNIQUE
    #         if seq_id in target_dict[gene_id]:
    #             logging.error("Target sequence is already in the list for '{0:s}' gene".format(gene_id))
    #             raise EXONtoolsError("Reference parsing error")
    #     except KeyError:
    #         target_dict[gene_id] = {}       # create empty dictionary (with GENE KEY) for selected gene

    #     # ADD ALL USEFUL INFO ABOUT THIS PARTICULAR SEQUENCE
    #     target_dict[gene_id].update({seq_id: {"id": seq_id, "seq": refseq.seq, "length": len(refseq.seq), "name": ref_name, "organism": org_name, "description": seq_info, "gene": gene_id}})

    # # CHOOSING THE LONGEST SEQUENCE FOR EACH GENE IF MULTIPLE WERE ASSIGNED
    # for i, gene in enumerate(sorted(target_dict.keys(), key=natural_sort)):
    #     prot_dict = sorted(target_dict[gene], key=lambda x: target_dict[gene][x]['length'], reverse=True)[0]
    #     target_dict["target" + str(i + 1)] = target_dict[gene][prot_dict]
    #     del target_dict[gene]

    logging.info("Total {0:d} reference sequences were parsed".format(len(target_dict)))

    targref.add_dict(target_dict)
    targref.add_path(reference)


def plustominus(sign):
    if sign == "+":
        return "-"
    elif sign == "-":
        return "+"
    else:
        return "."


def annotator_pars(query, target, evalue, cluster, minlen, similarity, overlap, database, custom, orf, translating, gencode, oneway, norename, organism):
    """Testing annotator parameters"""

    typedict = {"prot": "peptide", "nucl": "nucleotide"}
    gname = loadcodes(gencode)['name']

    if orf and translating:
        logging.error("'--orf' and '--translate' options cannot be used together in one EXONtools command")
        raise EXONtoolsError("Argument error")

    if orf and query == 'prot':
        logging.error("'--orf' option cannot be used with peptide query sequences")
        raise EXONtoolsError("Argument error")

    if translating and query == 'prot':
        logging.error("'--translate' option cannot be used with peptide query sequences")
        raise EXONtoolsError("Argument error")

    logging.debug("Testing annotator parameters")
    logging.info("Query contigs are defined as {0:s} sequences ('{1:s}' type)".format(typedict[query], query))
    logging.info("References are defined as {0:s} sequences ('{1:s}' type)".format(typedict[target], target))
    logging.info("NCBI genetic code --> {0:d} ('{1:s}')".format(gencode, gname))

    if query == 'nucl' and target == "nucl" and orf:
        logging.warning("All query contigs will be trimmed to long ORFs sequences")
        logging.info("The current annotation analysis will be performed using 'BLASTN' alignment")
    elif query == 'nucl' and target == "nucl" and translating:
        logging.warning("All query contigs will be trimmed to long ORFs and translated to AA sequences")
        logging.info("The current annotation analysis will be performed using 'TBLASTN' alignment")
    if query == 'nucl' and target == "prot" and orf:
        logging.warning("All query contigs will be trimmed to long ORFs sequences")
        logging.info("The current annotation analysis will be performed using 'BLASTX' alignment, frame 1")
    elif query == 'nucl' and target == "prot" and translating:
        logging.warning("All query contigs will be trimmed to long ORFs and translated to AA sequences")
        logging.info("The current annotation analysis will be performed using 'BLASTP' alignment")
    elif query == 'nucl' and target == 'prot':
        logging.info("The current annotation analysis will be performed using 'BLASTX' alignment")
    elif query == 'prot' and target == 'prot':
        logging.info("The current annotation analysis will be performed using 'BLASTP' alignment")
    elif query == 'prot' and target == 'nucl':
        logging.info("The current annotation analysis will be performed using 'TBLASTN' alignment")
    else:
        logging.info("All query contigs will be annotated using 'BLASTN' alignment")

    if evalue < 0:
        logging.error("The e-value argument cannot be negative")
        raise EXONtoolsError("Argument value error")
    else:
        logging.info("E-value alignment threshold is set to '{0:s}'".format(str(evalue)))

    if cluster < 0.8 or cluster > 1:
        logging.error("Sequence identity for contig filtration cannot be lower than 0.8 or greater than 1")
        raise EXONtoolsError("Argument value error")

    if minlen < 0:
        logging.error("Minumum contig length cannot have a negative value")
        raise EXONtoolsError("Argument value error")

    if similarity < 0 or similarity > 1:
        logging.error("Alignment similarity value cannot be lower than zero or greater than 1 ")
        raise EXONtoolsError("Argument value error")

    if overlap < 0:
        logging.error("Alignment length threshold cannot be less than 0")
        raise EXONtoolsError("Argument value error")

    if database == "custom" and custom:
        logging.warning("The reference database will be parsed according to user-defined custom format")
    elif database == "custom" and not custom:
        logging.error("The reference database is chosen to be in custom format but no custom format was provided")
        logging.error("Use '--grep' command option to set your custom format for parsing the database")
        raise EXONtoolsError("Argument value error")
    else:
        logging.info("The format of the provided reference file(s) was defined as '{0:s}'".format(database.upper()))

    if norename:
        logging.warning("All annotated contigs will retain their original names")
    else:
        logging.warning("All annotated contigs will be renamed according to their annotations")

    if oneway:
        logging.warning("All negative strand contigs will be automatically converted to positive strand")

    if organism:
        logging.warning("'{0:s}' will be used as organism identifier in all annotations".format(organism))

    logging.debug("Annotation parameter test passed: OK")
