# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.

import logging
import logging.config
import json
import os
import sys

# OUTPUT FORMAT SETTINGS
LOGFORMAT = '%(asctime)s EXONtools %(levelname)s -> %(module)s.py: %(message)s'
DATEFMT = "%m/%d/%y %H:%M:%S"


def setup_logger(level="INFO", logmode=False, quietmode=False):
    """This function initialize the logger for each EXONtools module."""

    loglevel = getattr(logging, level)
    if logmode:
        logpath = "EXONtools.log"
    else:
        logpath = False

    logger = logging.getLogger(__name__)
    logging.getLogger
    logger.setLevel(loglevel)
    logformat = logging.Formatter(LOGFORMAT, DATEFMT)

    if logmode:
        file_handler = logging.FileHandler(logpath)
        file_handler.setFormatter(logformat)
        logger.addHandler(file_handler)

    if quietmode:
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setLevel(logging.ERROR)
        stream_handler.setFormatter(logformat)
        logger.addHandler(stream_handler)
    else:
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setLevel(loglevel)
        stream_handler.setFormatter(logformat)
        logger.addHandler(stream_handler)

    return logger


def new_logger_settings(worker=None):
    """The function gets logging settings from json and runs new logger"""

    # jsonpath = os.path.join(os.path.split(sys.argv[0])[0], "jsons", "logging.json")
    jsonpath = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "jsons", "logging.json")

    with open(jsonpath, 'r') as logconfigfile:
        data = json.load(logconfigfile)
    if worker:
        data['logging']['formatters']['multi'] = {}
        data['logging']['formatters']['multi']['format'] = data['logging']['formatters']['simple']['format'].replace('->', '-> Worker-' + str(worker) + ' in')
        data['logging']['formatters']['multi']['datefmt'] = data['logging']['formatters']['simple']['datefmt']
        data['logging']['formatters']['multi']['class'] = data['logging']['formatters']['simple']['class']
        for handler in data['logging']['handlers']:
            data['logging']['handlers'][handler]['formatter'] = 'multi'
    return logging.config.dictConfig(data['logging'])


def log2json(level="INFO", logmode=False, quietmode=False, disable_existing_loggers=False):
    """Writes all logging parameters to json file"""

    # jsonpath = os.path.join(os.path.split(sys.argv[0])[0], "jsons", "logging.json")
    jsonpath = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "jsons", "logging.json")

    logging_dict = {
        "logging": {
            "version": 1,
            "disable_existing_loggers": disable_existing_loggers,
            "level": level,
            "formatters": {
                "simple": {
                    "class": "logging.Formatter",
                    "datefmt": DATEFMT,
                    "format": LOGFORMAT
                },
            },
            "handlers": {
                "console": {
                    "level": level,
                    "class": "logging.StreamHandler",
                    "formatter": "simple",
                    "stream": "ext://sys.stdout",
                }
            },
            "loggers": {},
            "root": {
                "handlers": ["console"],
                "level": "DEBUG"
            }
        }
    }

    if logmode:
        logging_dict["logging"]["handlers"].update({
            "file_handler": {
                "level": level,
                "class": "logging.handlers.WatchedFileHandler",
                "formatter": "simple",
                "filename": "EXONtools.log",
                "mode": "a",
                "encoding": "utf-8"
            }
        })
        logging_dict["logging"]['root']["handlers"].append("file_handler")

    if quietmode:
        logging_dict["logging"]["handlers"]["console"]["level"] = "ERROR"

    with open(jsonpath, 'w') as logconfigfile:
        json.dump(logging_dict, logconfigfile, indent=2)


def args2json(args):
    jsonpath = os.path.join(os.path.split(sys.argv[0])[0], "jsons", "lastcompars.json")

    arg_dict = {x: args.__dict__[x] for x in args.__dict__ if x != "command"}

    com_name = arg_dict.pop("action")
    jsondict = {com_name: arg_dict}
    with open(jsonpath, 'w') as argconfigfile:
        json.dump(jsondict, argconfigfile, indent=2)
