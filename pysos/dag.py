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
import copy

from .utils import env

#
# DAG design:
#
# 1. Nodes are jobs to be completed currently they are the
#    same as steps.
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
    def __init__(self, node_name, input_targets=[], depends_targets=[], output_targets=[]):
        self._node_id = node_name
        self._input_targets = [] if input_targets is None else copy.copy(input_targets)
        self._depends_targets = [] if depends_targets is None else copy.copy(depends_targets)
        self._output_targets = [] if output_targets is None else copy.copy(output_targets)
        #env.logger.error('Note {}: Input: {} Depends: {} Output: {}'.format(self._node_id, self._input_targets,
        #      self._depends_targets,  self._output_targets)) 

    def __repr__(self):
        return self._node_id

    def depends_on(self, node):
        return any(x in node._output_targets for x in self._input_targets) or \
            any(x in node._output_targets for x in self._depends_targets) 

class SoS_DAG(nx.DiGraph):
    def __init__(self):
        nx.DiGraph.__init__(self)

    def add_step(self, node_name, input_targets, depends_targets, output_targets):
        self.add_node(SoS_Node(node_name, input_targets, depends_targets, output_targets))

    def connect(self):
        '''Connect nodes according to status of targets'''
        # right now we do not worry about status of nodes
        #
        # connecting the output to the input of other nodes
        #
        # NOTE: This is implemented in the least efficient way just for
        # testing. It has to be re-implemented.
        # 
        # refer to http://stackoverflow.com/questions/33494376/networkx-add-edges-to-graph-from-node-attributes
        # 
        # for some code using attributes
        for node_i in self.nodes():
            for node_j in self.nodes():
                if node_i == node_j:
                    continue
                if node_i.depends_on(node_j):
                    self.add_edge(node_j, node_i)
                if node_j.depends_on(node_i):
                    self.add_edge(node_i, node_j)

    def write_dot(self, filename):
        try:
            wd = nx.drawing.nx_pydot.write_dot
            wd(self, filename)
        except Exception as e:
            env.logger.error('Failed to call write_dot: {}'.format(e))
