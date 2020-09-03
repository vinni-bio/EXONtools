# ENCODING: UTF-8


# This file was created by Kirill Vinnikov on August 10, 2020
# Copyright 2020 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root directory of the EXONtools package.

# The general idea of assemblation approach and the cooresponding
# perl script was provided in Bi, K., Vanderpool, D., Singhal, S.,
# Linderoth, T., Moritz, C. and Good, J.M., 2012. Transcriptome-based
# exon capture enables highly cost-effective comparative genomic data
# collection at moderate evolutionary scales. BMC Genomics, 13: 403.

from __future__ import print_function, division
import os
import re
import logging
import csv
import pdb

from mains.EXT_prog import EXTprogram
from mains.EXT_IO import getinput, output, makenewdir
from mains.EXT_executor import executor
from mains.EXT_worker import worker
from mains.EXT_parallel import hard_worker, create_pool, close_pool, run_instance, run_executor, set_threads
from mains.EXT_errors import EXONtoolsError
from utils.sorting import natural_sort
from utils.seqIO import SeqIO


class assemblator(EXTprogram):
    """This program merges assemblies and reduces redundancy of the final assemblies"""

    name = "assemblator"

    def execute_program(self):
        args = self.args
        self.reassemble(args.inpath, args.outdir, args.cluster, args.parsing, args.overlap, args.decrement, args.minlen, args.repeats, args.similarity)

    def reassemble(self, inpath, outdir, clusterid, parsing, overlap, decrement, minlen, repeats, similarity):

        if assemblator.debug:
            pdb.set_trace()

        executor.setconfig("blat", "cap3", "cd-hit-est")

        if assemblator.debug:
            pdb.set_trace()

        clusterid = round(clusterid, 2)
        similarity = round(similarity, 2)
        assemblator_pars(clusterid, parsing, overlap, minlen, similarity, decrement, repeats)

        if assemblator.debug:
            pdb.set_trace()

        # SET DRY RUN AND DEBUGGING MODES FOR SUBCLASSES
        assemblator.run_dry(SeqIO, getinput, output, makenewdir, executor)
        assemblator.set_debug(executor)

        # TESTING INPUT PARAMETER VALUES
        mem = int(assemblator.memory * 1000)

        # GET INPUT FILES
        getinput.format(['.fa', '.fasta'])
        FileList = getinput(inpath).files

        # MAKE OUTPUT DIRECTORY
        output(outdir)

        # MAKE TMP DIRECTORY
        tmpdir = makenewdir(name="tmp", fullname="temporary")

        # MAKE LOG DIR
        log_dir = makenewdir(name="dependency_logs", fullname="LOG")
        logging.info("Find all dependency program logs in the '{0:s}' folder".format(os.path.basename(log_dir.path)))

        logging.debug("IO settings: OK")

        logging.debug("Ready to start assemblation analysis")

        if assemblator.debug:
            pdb.set_trace()

        # CREATING WORKING VARS
        libs = {}
        working_paths = {}
        Nstage = 0

        # 0. MERGING
        logging.info("**********************************************************************")
        logging.info("STAGE {0:d}. MERGING ASSEMBLIES".format(Nstage))
        logging.info("**********************************************************************")
        Nstage += 1
        concat = makenewdir(name=os.path.join(tmpdir.path, "concat"), fullname='concatenated assemblies')
        working_paths = {}
        for f in FileList:
            if parsing:
                libname = os.path.basename(f).split("_")[0].split(".")[0]
            else:
                libname = "CONSENSUS"
            if libname not in libs:
                libs[libname] = {"contigs": 0}
            outpath = os.path.join(concat.path, libname + "_conc.fa")
            working_paths[libname] = outpath
            mergefile = SeqIO(f)

            for x in mergefile.read():
                libs[libname]["contigs"] += 1
                with open(outpath, 'a') as outfile:
                    outfile.write(">contig{0:d}\n{1:s}\n".format(libs[libname]["contigs"], x.seq))

        logging.info("{0:d} libraries in total will be reassembled with 'assemblator' algorithm".format(len(libs)))
        for lib in sorted(libs.keys(), key=natural_sort):
            logging.debug("{0:s} with {1:d} merged contigs".format(lib, libs[lib]["contigs"]))
        logging.debug("Merging assemblies is finished: OK")

        if assemblator.debug:
            pdb.set_trace()

        # 1. CLUSTERTING (100% IDENTITY)
        logging.info("**********************************************************************")
        logging.info("STAGE {0:d}. PREPOCESSING".format(Nstage))
        Nstage += 1
        logging.info("**********************************************************************")
        logging.info("Filtering contig duplicates (100% identity) and contigs shorter than {0:d}bp".format(minlen))
        clust100 = makenewdir(name=os.path.join(tmpdir.path, "clust100"), fullname='temporary cluster')
        for libname in sorted(libs.keys(), key=natural_sort):
            outpath = os.path.join(clust100.path, libname + "_clust100.fa")
            outlog = os.path.join(log_dir.path, libname + "_preprocessing.log")
            outclust = outpath + ".clstr"
            logging.debug("Clustering contigs in '{0:s}' library with 'cd-hit-est'".format(libname))
            run_executor(executor(
                program="cd-hit-est",
                params=[working_paths[libname], outpath, 1.00, minlen, mem, assemblator.threads],
                conditions={"positive": [assemblator.threads], "pathexists": [working_paths[libname]]},
                custom_arg_string=assemblator.extra + " > " + outlog + " 2>&1")
            )
            working_paths[libname] = outpath
            clusters = parse_cluster(outclust)
            libs[libname].update({"clust100": len(clusters)})
            logging.info("Total {0:d} contigs retained in '{1:s}' library".format(len(clusters), libname))
        concat.delete()
        logging.debug("Clustering with 100% identity is finished: OK")

        if assemblator.debug:
            pdb.set_trace()

        # 2. CLUSTERING AND ASSEMBLING CONTIGS
        if clusterid < 1:
            logging.info("**********************************************************************")
            logging.info("STAGE {0:d}. INITIAL CLUSTERIZATION".format(Nstage))
            Nstage += 1
            clustdir = makenewdir(name=os.path.join(tmpdir.path, "clustdir"), fullname='temporary cluster')

            logging.debug("Grouping and assembling contigs within each library using CDHIT and CAP3")
            for libname in sorted(libs.keys(), key=natural_sort):
                count = 0
                nclust = 0
                cap3dir = makenewdir(os.path.join(clustdir.path, "tmpass"))
                outpath = os.path.join(clustdir.path, libname + "_cluster.fa")
                outlog = os.path.join(log_dir.path, libname + "_cluster.log")
                outclust = outpath + ".clstr"

                logging.info("**********************************************************************")
                logging.info("Grouping contigs with {0:d}% identity in '{1:s}' library".format(int(clusterid * 100), libname))
                run_executor(executor(
                    program="cd-hit-est",
                    params=[working_paths[libname], outpath, clusterid, minlen, mem, assemblator.threads],
                    conditions={"positive": [assemblator.threads], "pathexists": [working_paths[libname]]},
                    custom_arg_string=assemblator.extra + " > " + outlog + " 2>&1")
                )

                try:
                    os.remove(outpath)
                except IOError:
                    logging.error("Clusterization error. No output file was produced for '{0:s}' library".format(libname))
                    raise EXONtoolsError("Clusterization output error")

                logging.info("Reading resulted clusters within '{0:s}' library".format(libname))
                clusters = parse_cluster(outclust)
                contigs = convert_clusters(clusters)

                conts100 = SeqIO(working_paths[libname])
                outpath = os.path.join(clustdir.path, libname + "_clustered.fa")
                cap3files = []
                with open(outpath, 'a') as outfile:
                    for contig in conts100.read():
                        if len(clusters[contigs[contig.name]]) > 1:
                            nclust += 1
                            inpathcap = os.path.join(cap3dir.path, libname + "_" + str(contigs[contig.name]) + ".fa")
                            cap3files.append(inpathcap)
                            with open(inpathcap, "a") as capfile:
                                capfile.write(">{0:s}\n{1:s}\n".format(contig.name, contig.seq))
                        else:
                            count += 1
                            outfile.write(">contig{0:d}\n{1:s}\n".format(count, contig.seq))

                del clusters
                del contigs
                working_paths[libname] = outpath
                cap3files = list(set(cap3files))

                logging.info("{2:d} contigs from '{1:s}' library are grouped into {3:d} clusters with {0:d}% identity level".format(int(clusterid * 100), libname, nclust, len(cap3files)))

                if cap3files:
                    logging.info("Assembling {0:d} contigs in '{1:s}' library using 'CAP3' program".format(nclust, libname))
                    TASKS = []
                    for f in cap3files:
                        TASKS.append(worker(cap3_processor, [f]))

                    processes_requested = set_threads("CAP3 processor", len(TASKS), assemblator.threads)
                    pool = create_pool(processes_requested)
                    jobs = hard_worker(run_instance, TASKS, pool)
                    close_pool(pool)
                    with open(outpath, 'a') as outfile:
                        for result in jobs:
                            for x in result:
                                count += 1
                                outfile.write(">contig{0:d}\n{1:s}\n".format(count, x))
                    del jobs
                    del TASKS

                if not os.listdir(cap3dir.path):
                    cap3dir.delete()
                else:
                    raise EXONtoolsError("CAP3 foder is not empty!")
                libs[libname].update({"clusters": count})
                logging.info("Total {0:d} contigs retained in '{1:s}' library after clusterization step".format(count, libname))
        logging.debug("Assembling of contig clusters with {0:d}% identity finished: OK".format(int(clusterid * 100)))

        if assemblator.debug:
            pdb.set_trace()

        # 3. REPETATIVE SELF-ALIGNMENT (BLAT)
        logging.info("**********************************************************************")
        logging.info("STAGE{0:d}. CONSENSUS GENERATION".format(Nstage))
        Nstage += 1
        logging.info("**********************************************************************")
        logging.info("BLAT clusterization will make {0:d} iterations on each library".format(repeats))
        logging.info("-----starting contig similarity threshold = {0:0.2f}".format(similarity))
        logging.info("-----alignment/overlap threshold = {0:0.2f}".format(overlap))
        logging.info("-----similarity decrement within each cycle = {0:0.2f}".format(decrement))

        for libname in sorted(libs.keys(), key=natural_sort):

            file_ids = []
            for step in list(range(repeats)):
                blatinput = working_paths[libname]
                logging.info("**********************************************************************")
                logging.info("BLAT clusterization of '{1:s}' library: cycle #{0:d}".format(step + 1, libname))
                outlog = os.path.join(log_dir.path, libname + "_BLAT_STEP" + str(step + 1) + ".log")
                if step == 0:
                    identity = similarity
                elif 0 < step <= repeats:
                    identity = round(identity - decrement, 2)
                else:
                    logging.error("Wrong BLAT clusterization settings")
                    raise EXONtoolsError("Assemblator error in BLAT procedure")

                filename = "_blat" + str(int(identity * 100))
                file_ids.append(filename)

                blatpath = os.path.join(tmpdir.path, filename[1:])

                if not os.path.exists(blatpath):
                    makenewdir(blatpath)

                logging.info("Contig similarity is set to: {0:0.2f}".format(identity))

                outblatpath = os.path.join(blatpath, libname + filename + ".out")

                run_executor(executor(
                    program="BLAT",
                    params=[blatinput, blatinput, assemblator.threads, outblatpath],
                    conditions={"positive": [assemblator.threads], "pathexists": [blatinput]},
                    custom_arg_string=" >> " + outlog + " 2>&1")
                )

                clusters = parse_blat(outblatpath, overlap, identity)
                contigs = convert_clusters(clusters)

                blatfile = SeqIO(blatinput)
                outpath = os.path.join(blatpath, libname + filename + ".fa")
                cap3dir = makenewdir(os.path.join(blatpath, "tmpass"))
                cap3files = []
                count = 0
                nclust = 0
                with open(outpath, 'a') as outfile:
                    for contig in blatfile.read():
                        try:
                            inpathcap = os.path.join(cap3dir.path, libname + "_cluster" + str(contigs[contig.name]) + ".fa")
                            nclust += 1
                            cap3files.append(inpathcap)
                            with open(inpathcap, "a") as capfile:
                                capfile.write(">{0:s}\n{1:s}\n".format(contig.name, contig.seq))
                        except KeyError:
                            count += 1
                            outfile.write(">contig{0:d}\n{1:s}\n".format(count, contig.seq))

                cap3files = list(set(cap3files))
                logging.info("{2:d} contigs from '{1:s}' library are grouped into {3:d} clusters with {0:d}% identity level".format(int(identity * 100), libname, nclust, len(cap3files)))

                del clusters
                del contigs
                working_paths[libname] = outpath

                if cap3files:
                    TASKS = []
                    for f in cap3files:
                        TASKS.append(worker(cap3_processor, [f]))

                    processes_requested = set_threads("CAP3 processor", len(TASKS), assemblator.threads)
                    pool = create_pool(processes_requested)
                    jobs = hard_worker(run_instance, TASKS, pool)
                    close_pool(pool)
                    with open(outpath, 'a') as outfile:
                        for result in jobs:
                            for x in result:
                                count += 1
                                outfile.write(">contig{0:d}\n{1:s}\n".format(count, x))
                    del jobs
                    del TASKS

                if not os.listdir(cap3dir.path):
                    cap3dir.delete()
                else:
                    raise EXONtoolsError("CAP3 foder is not empty!")

                libs[libname].update({filename[1:]: count})

                logging.info("Total {0:d} contigs retained in '{1:s}' library after #{2:d} cycle".format(count, libname, step + 1))
            logging.info("**********************************************************************")

            logging.debug("BLAT clusterization step {0:d} finished: OK".format(step + 1))

        logging.debug("FINISHED ALL BLAT ALIGNMENTS: OK")

        if assemblator.debug:
            pdb.set_trace()

        # 4. FINAL CD-HIT-EST TO FILTER 99% SIMILAR CONTIGS
        logging.info("**********************************************************************")
        logging.info("STAGE{0:d}. FINAL CLUSTERIZATION".format(Nstage))
        Nstage += 1
        for libname in sorted(libs.keys(), key=natural_sort):
            logging.info("**********************************************************************")
            logging.info("Performing the final clusterization with {0:d}% identity in '{1:s}' library".format(int(100 * clusterid), libname))
            path_final = os.path.join(output.path, libname + assemblator.suffix + ".fa")
            outlog = os.path.join(log_dir.path, libname + "_clust99_final.log")
            run_executor(executor(
                program="cd-hit-est",
                params=[working_paths[libname], path_final, clusterid, minlen, mem, assemblator.threads],
                conditions={"positive": [assemblator.threads], "pathexists": [working_paths[libname]]},
                custom_arg_string=assemblator.extra + " > " + outlog + " 2>&1")
            )
            clusters = parse_cluster(path_final + '.clstr')
            libs[libname].update({"FINAL": len(clusters)})
            logging.info("Total {0:d} contigs were found in the final assembly '{1:s}'".format(len(clusters), os.path.basename(path_final)))
        logging.debug("Final clusterization finished: OK")
        logging.info("**********************************************************************")
        # identity = similarity

        if assemblator.debug:
            pdb.set_trace()

        if assemblator.stats:
            statsfile = os.path.join(output.path, "consensus_assembly_stats.csv")
            logging.info("All summary data on contig processing in all assembling will be saved to:".format(os.path.basename(statsfile)))
            header = ["RRL", "MERGED", "UNIQUE", "ASSEMBLED"]
            for h in file_ids:
                header.append(h[1:].upper())
            header.append("FINAL_CONTIGS")

            if not assemblator.dryrun:
                with open(statsfile, 'w') as statout:
                    csv_writer = csv.writer(statout)
                    csv_writer.writerow(header)
                    for lib in sorted(libs.keys(), key=natural_sort):
                        csvrow = [lib, libs[lib]['contigs'], libs[lib]['clust100'], libs[lib]['clusters']]
                        for x in file_ids:
                            csvrow.append(libs[lib][x[1:]])
                        csvrow.append(libs[lib]["FINAL"])
                        csv_writer.writerow(csvrow)

        if not assemblator.keeptmp:
            tmpdir.delete()


def parse_cluster(inpath):
    """This function collects contig IDs within each cluster and returns the dictionary"""

    clusters = {}
    if not assemblator.dryrun:
        try:
            with open(inpath, 'r') as infile:
                for line in infile:
                    line = line.strip()
                    if line.startswith(">"):
                        mask = re.search("^>(Clus.*)$", line)
                        cluster = mask.group(1)
                        cluster = cluster.replace(" ", "")
                        clusters[cluster] = []
                    else:
                        mask = re.search(">(contig\d+)", line)
                        if mask:
                            clusters[cluster].append(mask.group(1))
        except IOError:
            logging.error("Cluster file {0:s} does not exist".format(inpath))
            raise EXONtoolsError("Non-existing path error")
        os.remove(inpath)
    return clusters


def convert_clusters(clusdict):
    """Convert values to clusters and keys to contignames"""

    newdict = {}
    count = 0

    for clustid, conts in clusdict.items():
        for contname in conts:
            count += 1
            newdict[contname] = clustid

    if count != len(newdict):
        raise EXONtoolsError("Convert clusters error")

    return newdict


def readvalidseq(inpath):
    with open(inpath, 'r') as infile:
        seqlist = []
        line = infile.readline().strip()
        while(line):
            if line.startswith(">"):
                seqline = ""
                line = infile.readline().strip()
                while(line):
                    seqline = seqline + line
                    line = infile.readline().strip()
                    if line is None or line.startswith(">"):
                        break
                seqlist.append(seqline)
            else:
                line = infile.readline().strip()
        return seqlist


def writevalidseq(outpath, SEQS):
    with open(outpath, 'w') as outfile:
        for i, x in enumerate(SEQS):
            outfile.write(">contig{0:d}\n{1:s}\n".format(i + 1, x))


def cap3_processor(inpath, maxseqs=700):
    """This function runs CAP3 and parses the output"""

    inseqs = readvalidseq(inpath)
    seqnum = len(inseqs)
    if seqnum > maxseqs:
        singlets = []
        contigs = []
        subsets = [inseqs[i * maxseqs:(i + 1) * maxseqs] for i in range((seqnum + maxseqs - 1) // maxseqs)]
        count = 0
        for seqs in subsets:
            count += 1
            outpath = inpath + ".part" + str(count) + ".fasta"
            writevalidseq(outpath, seqs)
            run_executor(executor(
                program="CAP3",
                params=[outpath],
                conditions={},
                custom_arg_string=" > /dev/null 2>&1",
                quiet=True)
            )
            contigs = contigs + readvalidseq(outpath + ".cap.contigs")
            singlets = singlets + readvalidseq(outpath + ".cap.singlets")
        del subsets

        if singlets:
            outpath = inpath + "_singlets.fasta"
            writevalidseq(outpath, singlets)
            run_executor(executor(
                program="CAP3",
                params=[outpath],
                conditions={},
                custom_arg_string=" > /dev/null 2>&1",
                quiet=True)
            )
            singlets = readvalidseq(outpath + ".cap.singlets")

        if contigs:
            outpath = inpath + ".contigs.fasta"
            writevalidseq(outpath, contigs)
            run_executor(executor(
                program="CAP3",
                params=[outpath],
                conditions={},
                custom_arg_string=" > /dev/null 2>&1",
                quiet=True)
            )
            contigs = readvalidseq(outpath + ".cap.contigs")

        finalseqs = contigs + singlets
        writevalidseq(inpath, finalseqs)
        del contigs
        del singlets
        del finalseqs

    run_executor(executor(
        program="CAP3",
        params=[inpath],
        conditions={},
        custom_arg_string=" > /dev/null 2>&1",
        quiet=True)
    )

    SEQS = readvalidseq(inpath + ".cap.contigs") + readvalidseq(inpath + ".cap.singlets")

    pardir = os.path.dirname(os.path.realpath(inpath))
    pattern = os.path.basename(inpath)
    for x in os.listdir(pardir):
        if x.startswith(pattern):
            os.remove(os.path.join(pardir, x))
    return SEQS


def parse_blat(inpath, threshold, similarity):
    """THIS FUNCTION SEARCHES FOR CONTIG CLUSTERS IN THE BLAT ALIGNMENT"""

    clusters = {}                                               # dictionary {clusters : values}
    revClusters = {}                                            # dictionary {contigs : values}
    count = 0                                                   # cluster counter

    if not assemblator.dryrun:
        try:
            with open(inpath, 'r') as infile:
                for line in infile:                                     # read BLAT output in psl format
                    psl = pslline(line.strip())
                    if not psl.selfcheck() and psl.alfrac() >= threshold and psl.simfrac() >= similarity:
                        query_id = None                              # create blank query variable for cluster ID
                        target_id = None                             # create blank target variable for cluster ID
                        if psl.query in revClusters:                 # if query in revClusters
                            query_id = revClusters[psl.query]        # assign cluster ID to query_id
                        if psl.target in revClusters:                # if target in revClusters
                            target_id = revClusters[psl.target]      # assign cluster ID to target_id
                        if query_id and target_id:                   # if both were found in revClusters
                            if query_id != target_id:                # as different cluster IDs
                                for i in clusters[target_id]:        # then transfer all contigs from target cluster
                                    clusters[query_id].append(i)     # to query cluster
                                    revClusters[i] = query_id        # change cluster IDs for these contigs in revCluster
                                del clusters[target_id]              # and delete the target cluster
                        elif query_id:                               # if only query in revClusters
                            clusters[query_id].append(psl.target)    # append target contig to corresponding cluster
                            revClusters[psl.target] = query_id       # add target contig to revCluster
                        elif target_id:                              # if only target in revClusters
                            clusters[target_id].append(psl.query)    # append query contig to corresponding cluster
                            revClusters[psl.query] = target_id       # add query contig to revCluster
                        else:
                            count += 1                               # count clusters
                            clusters[count] = []                     # creates empty cluster with count ID
                            clusters[count].append(psl.query)        # appends query contig to the new cluster
                            clusters[count].append(psl.target)       # appends target contig to the new cluster
                            revClusters[psl.query] = count           # adds query contig to revClusters with cluster ID
                            revClusters[psl.target] = count          # adds target contig to revClusters with cluster ID
        except IOError:
            logging.error("The BLAT output file cannot be opened from the provided path")
            raise EXONtoolsError("Non-existing path error in BLAT parse function")
        os.remove(inpath)
    return clusters


class pslline(object):
    """Parses the psl line"""

    def __init__(self, line):

        if isinstance(line, str):
            pslist = line.split("\t")
            pslist = [int(x) if x.isdigit() else x for x in pslist]

            self.match = pslist[0]                   # Number of matching bases that aren't repeats.
            self.mismatch = pslist[1]                # Number of bases that don't match.
            # self.nbases = pslist[3]                  # Number of 'N' bases.
            # self.qnins = pslist[4]                   # Number of inserts in query.
            # self.qbins = pslist[5]                   # Number of bases inserted into query.
            # self.tnins = pslist[6]                   # Number of inserts in target.
            # self.tbins = pslist[7]                   # Number of bases inserted into target.
            # self.strand = pslist[8]                  # defined as + (forward) or - (reverse) for query strand.
            self.query = pslist[9]                   # Query sequence name.
            self.qsize = pslist[10]                  # Query sequence size.
            self.qstart = pslist[11]                 # Alignment start position in query.
            self.qend = pslist[12]                   # Alignment end position in query.
            self.target = pslist[13]                 # Target sequence name.
            self.tsize = pslist[14]                  # Target sequence size.
            self.tstart = pslist[15]                 # Alignment start position in target.
            self.tend = pslist[16]                   # Alignment end position in target.
            self.allen = self.match + self.mismatch
        else:
            logging.error("PSL line should be in string format")
            raise EXONtoolsError("PSL line format error")

    def selfcheck(self):
        """Checks self-alignents"""

        if self.query == self.target:
            return True
        else:
            return False

    def alfrac(self):
        """return fraction of alignment length in overlapping region"""

        if self.qstart >= self.tstart and self.qsize - self.qend >= self.tsize - self.tend:     # if target within query
            return self.allen / self.tsize
        elif self.qstart >= self.tstart and self.qsize - self.qend < self.tsize - self.tend:    # if target extends to the right from query
            return self.allen / (self.tend + self.qsize - self.qend)
        elif self.qstart < self.tstart and self.qsize - self.qend >= self.tsize - self.tend:    # if target extends to the left from query
            return self.allen / (self.qend + self.tsize - self.tend)
        elif self.qstart < self.tstart and self.qsize - self.qend < self.tsize - self.tend:     # if query within target
            return self.allen / self.qsize
        else:
            logging.error("Something went wrong in BLAT alignment overlap comparisons...")
            raise EXONtoolsError("Assemblator error in BLAT procedure")

    def simfrac(self):
        """return fraction of matched bases in alignment"""

        return self.match / self.allen


def assemblator_pars(cluster, parsing, overlap, minlen, similarity, decrement, repeats):
    """This function verifies all paramaters provided for assemblator analysis"""

    logging.debug("Checking assemblator parameters")
    if assemblator.extra:
        logging.warning("The following extra arguments will be added to CD-HIT-EST command line:")
        logging.warning(assemblator.extra)
    if repeats < 0:
        logging.error("The number of BLAT cycle repeats must be equal or above 0")
        raise EXONtoolsError("Argument value error")
    if parsing:
        logging.warning("Input files will be grouped by library IDs to produce separate assemblies")
    if decrement > 0.1 or decrement < 0:
        logging.error("The decrement value cannot be greater than 0.1 or lower than 0")
        raise EXONtoolsError("Argument value error")
    if similarity - repeats * decrement < 0.6:
        logging.error("The final similarity threshold cannot be lower than 60%%")
        raise EXONtoolsError("Argument value error")
    if minlen < 0:
        logging.error("The minumum contig length cannot be smaller than 0")
        raise EXONtoolsError("Argument value error")
    if similarity < 0 or similarity > 1:
        logging.error("Starting similarity cannot be greater than 1 or lower than 0")
        raise EXONtoolsError("Argument value error")
    if overlap > 1 or overlap <= 0:
        logging.error("Overlap threshold cannot be more than 1, equal or less than 0")
        raise EXONtoolsError("Argument value error")
    if cluster > 1 or cluster < 0.9:
        logging.error("Initital sequence identity level cannot be more than 1, equal or less than 0.9")
        raise EXONtoolsError("Argument value error")

    logging.debug("Assemblator parameter check: OK")
