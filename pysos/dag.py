#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
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

import networkx as nx

#
# DAG design:
#
# 1. Nodes are jobs to be completed, including
#
#    1.a. Steps.
#    1.b. Tasks within steps.
#
#    The completion of steps depends on the completion of tasks
#    so we should see
#
#       task --> step
#
#    If there are multiple tasks within a step, the graph would look
#    like 
#
#       task 0 -> task 1 -> task 2 -> step
#   
#    if concurrent = False (default), or
#
#       task 0  ->
#       task 1  -> step
#       task 2  -> 
#
#    if concurrent = True.
#
# 2. Input, output and depedent files of nodes.
#
#    Each step has input, output and depdent files with their status.
#    The status of these files determines the status of ndoes.
#    There are several status for these targets
#
#    a. Target exists and without signature.
#       The files need to be processed and obtain its signature.
#
#    b. Target does not exists and with signature.
#       The files are assumed to exist for dependency check, but
#       are determined to be non-exist if they are needed for a 
#       task.
#
#    c. Target does not exists and without signature.
#       The target need to be generated and obtain its signature.
#
#    d. Target is undetermined. A node with undetermined input
#       will depend on all its previous steps.
#
#  3. Status of node.
#
#    a. A node does not depend on any input file.
#       It can be executed immediately.
#
#    b. A node depends on undetermined input file.
#       It depend on all previous steps and can only be processed
#       when all its previous steps are implemented.
#
#    c. A node depends on a previous file.
#       It can only be executed when the previous step is executed.
#       This is the regular case.
#
#    d. A node depends on a file that does not exist.
#       This is problematic. An error will be raised.
#
# 4. Edges.
#
#    a. Edges determines dependencies between nodes. They are
#       generated depending on the status of the nodes.
#    
#    b. Edges will be updated with the update of nodes.
#
# 5. Dynamic DAG.
#
#    a. DAG will consists of nodes and edges derived from 
#       status of files.
#
#    b. Completion of a task or step will update the status of
#       files. All edges related to updated targets will be updated.
#
#    c. Addition of node is currently not considered, but should
#       be allowed.
#
class SoS_Node(object):
    def __init__(self, node_id, input_targets=[], output_targets=[]):
        self._node_id = node_id
        self._input_targets = input_targets
        self._output_targets = output_targets

class SoS_DAG(nx.DiGraph):
    def __init__(self):
        nx.DiGraph.__init__(self)

    def add_step(self, node):
        self.add_node(node)

    def add_task(self, step, task):
        self.add_node(task)

