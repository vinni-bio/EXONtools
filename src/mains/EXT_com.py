# ENCODING: UTF-8

# This file was created by Kirill Vinnikov on August 10, 2018
# Copyright 2018 by Kirill Vinnikov. All rights reserved.

# This code is a part of the EXONtools distribution and governed
# by its license. Please see the LICENSE.txt file that should
# have been included in the root folder of the EXONtools package.


from __future__ import print_function, unicode_literals
import logging


class EXTcommand:
    """Represents a portable, executable Command.  EXTcommand have a list of supported_programs, a
        default_program, and an args object (object with parameters as attributes).  EXONcommands read the
        program attribute from the args object, and call the appropriate program (if supported), falling back to the
        default if no match is found.
    """
    supported_programs = []
    default_program = None
    command_name = ""

    def __init__(self, args):
        self.args = args

    def get_program(self, program):
        """Given a program name, searches for that Program in the list of supported_programs, and returns an instance
        of that Program initalized with the parameter object.  If not found, returns an instance of the default_program
        initalized with the parameter object.
        """
        for prog in self.supported_programs:
            if prog.name == program:
                logging.debug("The program '{}' was selected for running the analysis".format(prog.name))
                return prog(self.args)
        return self.default_program()

    def execute_command(self):
        """Executes the command.

        :return: Results of execution.
        """
        prog = self.get_program(self.args.program)
        logging.info("EXONtools is executing the command '{}' using the program '{}'".format(self.command_name, prog.name))

        return prog.execute_program()
