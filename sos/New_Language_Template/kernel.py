#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS for more information.
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#


###############################################################################
#                                                                             #
# To define new language for your particular language                         #
#                                                                             #
# 1. Define a class with                                                      #
#    * an __init__ function that accepts the SoS kernel.                      #
#    * an member kernel_name, which should be the kernel name for the         #
#      language that you can see from the output of command 'jupyter          #
#      kernerlspec list'                                                      #
#    * a background color in HTML hex format. This will be the background     #
#      color of the prompt area of the cells using this language. You can     #
#      set it to blank to stop highlighting the cells.                        #
#    * an init_statements that will be executed in the subkernel when the     #
#      kernel is started.                                                     #
#    * a function sos_to_lan that accepts an Python object and its name.      #
#      The function should return a new_name (can be the same if the name     #
#      is acceptable in LAN), and a statement to be executed in the subkernel #
#      to create the variable in the subkernel.                               #
#    * a function lan_to_sos that accepts a list of names from the subkernel. #
#      This function is responsible of getting a Python dictionary of name    #
#      (key) and object (value) for the specified items, and all variables    #
#      in the subkernel with names starting with sos. You will need to use    #
#      function get_response of the sos kernel to execute one or more         #
#      statements in the subkernel to achive this.                            #
# 2. Modify setup.py and add the language to the [sos_languages] section      #
#    of entry_points.                                                         #
#                                                                             #
# The following is a placeholder class for a new language. Please refer to    #
# implementation of existing languages and the 'Extending_SoS' section of     #
# the SoS documentation (vatlab.github.io/SOS/) for details.                  #
#                                                                             #
###############################################################################


# class sos_LAN:
#     def __init__(self, sos_kernel):
#         self.sos_kernel = sos_kernel
#         self.kernel_name = 'LAN'
#         self.background_color = '#XXXXXX'
#         self.init_statements = ''
#
#     def sos_to_lan(self, name, obj):
#         return name, ''
#
#     def lan_to_sos(self, items):
#         return {}



