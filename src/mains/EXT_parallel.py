# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

from __future__ import print_function, division
import sys
import signal
import multiprocessing as mp
import logging
from mains.EXT_logger import new_logger_settings
from mains.EXT_worker import worker
from time import sleep


def run_instance(worker_instance):
    """Task launcher. Runs a single worker instance"""
    return worker_instance.run()


def run_executor(program):
    """Task launcher. Runs a single 'Executor' instance."""
    return program.run_program()


def hard_worker(runner, TASKS, pool):
    """Runs all tasks and saves their results """
    try:
        results = pool.map_async(runner, TASKS).get(999999999)
        return results
    except KeyboardInterrupt:
        logging.exception("The current EXONtools run was interrupted by keyboard STOP command")
        kill_pool(pool)
    except AssertionError:
        logging.exception("EXONtools multiprocessing error")
        kill_pool(pool)
    except Exception:
        logging.error(sys.exc_info())
        logging.exception("EXONtools multiprocessing error")
        kill_pool(pool)


def create_pool(processes):
    """Creates pool for given number of processes"""
    pool_size = max(1, processes)
    pool = mp.Pool(pool_size, init_worker)
    logging.debug("Hiring {0:d} workers to do all the hard work: OK".format(pool_size))
    return pool


def close_pool(pool):
    """Cleans and closes the pool"""
    pool.close()
    pool.join()
    logging.debug("{0:d} jobs are done in total. Fire off all workers".format(worker.tasks_performed))
    worker.reset()


def kill_pool(pool):
    """Cleans and kills the pool with program exit"""
    logging.error("Encountered multiprocessing error! Closing the worker pool and exit...\n")
    logger = logging.getLogger('my-logger')
    logger.propagate = False
    pool.close()
    pool.terminate()
    pool.join()
    sys.exit()


def init_worker():
    """Handles signal from process interruption"""
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    worker_id = mp.current_process()._identity[0]
    sleep(worker_id / 2)
    new_logger_settings(worker=worker_id)


def set_threads(program, Ntasks, threads, userdef=False):
    """Returns the number threads and prints the output"""
    max_threads = mp.cpu_count()
    logging.debug("CPUs are found in your system: {0:d}".format(max_threads))
    logging.debug("CPUs requested by user: {0:d}".format(threads))
    if userdef:
        processes_requested = threads
    elif threads > max_threads:
        logging.info("The number of requested threads was changed to {0:d}".format(max_threads))
        threads = max_threads
    else:
        processes_requested = min(Ntasks, threads)
    logging.debug("{0:d} jobs are created in total".format(Ntasks))
    logging.debug("{1:d} CPU(s) will be used by '{0:s}'".format(program, processes_requested))
    return processes_requested
