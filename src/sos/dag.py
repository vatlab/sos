#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/vatlab/SOS for more information.
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
from collections import defaultdict
import copy
import pickle
import time
import fasteners

from .utils import env, ActivityNotifier, short_repr
from .eval import Undetermined
from .targets import file_target, sos_variable, textMD5, sos_step, sos_targets


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
    def __init__(self, step_uuid, node_name, node_index, input_targets=None, depends_targets=None,
        output_targets=None, context=None):
        self._step_uuid = step_uuid
        self._node_id = node_name
        self._node_index = node_index
        self._input_targets = Undetermined() if input_targets is None else copy.deepcopy(input_targets)
        self._depends_targets = [] if depends_targets is None else copy.deepcopy(depends_targets)
        self._output_targets = Undetermined() if output_targets is None else copy.deepcopy(output_targets)
        #env.logger.error('Note {}: Input: {} Depends: {} Output: {}'.format(self._node_id, self._input_targets,
        #      self._depends_targets,  self._output_targets))
        self._context = {} if context is None else copy.deepcopy(context)
        if '__completed__' not in self._context:
            self._context['__completed__'] = []
        self._status = None
        # unique ID to avoid add duplicate nodes ...
        self._node_uuid = textMD5(pickle.dumps((step_uuid, node_name, node_index,
                input_targets, depends_targets, output_targets,
                [] if context is None else [(k, sorted(list(context[k])) if isinstance(context[k], set) else context[k])
                    for k in sorted(context.keys())])))

    def __repr__(self):
        return self._node_id

    def show(self):
        print(
            f'{self._node_id} ({self._node_index}, {self._status}): input {self._input_targets}, depends {self._depends_targets}, output {self._output_targets}, context {self._context}')

class SoS_DAG(nx.DiGraph):
    def __init__(self, *args, **kwargs):
        nx.DiGraph.__init__(self, *args, **kwargs)
        # all_dependent files includes input and depends files
        self._all_dependent_files = defaultdict(list)
        self._all_output_files = defaultdict(list)

    def num_nodes(self):
        return nx.number_of_nodes(self)

    def add_step(self, step_uuid, node_name, node_index, input_targets, depends_targets,
        output_targets, context=None):
        node = SoS_Node(step_uuid, node_name, node_index, input_targets, depends_targets,
            output_targets, context)
        if node._node_uuid in [x._node_uuid for x in self.nodes()]:
            return
        # adding a step would add a sos_step target to met the depends on sos_step
        # requirement of some steps.
        self._all_output_files[sos_step(node_name.split(' ')[0])].append(node)
        if isinstance(input_targets, sos_targets) or (not isinstance(input_targets, Undetermined) and any(isinstance(x, sos_targets) for x in input_targets)):
            raise RuntimeError(f"input_targets should not be sos_targets {input_targets}")
        if isinstance(depends_targets, sos_targets) or (not isinstance(depends_targets, Undetermined) and any(isinstance(x, sos_targets) for x in depends_targets)):
            raise RuntimeError(f"depends_targets should not be sos_targets {depends_targets}")
        if isinstance(output_targets, sos_targets) or (not isinstance(output_targets, Undetermined) and any(isinstance(x, sos_targets) for x in output_targets)):
            raise RuntimeError(f"output_targets should not be sos_targets {output_targets}")
        if not isinstance(input_targets, (type(None), Undetermined)):
            for x in input_targets:
                if node not in self._all_dependent_files[x]:
                    self._all_dependent_files[x].append(node)
        if not isinstance(depends_targets, (type(None), Undetermined)):
            for x in depends_targets:
                if node not in self._all_dependent_files[x]:
                    self._all_dependent_files[x].append(node)
        if not isinstance(output_targets, (type(None), Undetermined)):
            for x in output_targets:
                if node not in self._all_output_files[x]:
                    self._all_output_files[x].append(node)
        if context is not None:
            for x in context['__changed_vars__']:
                if node not in self._all_output_files[sos_variable(x)]:
                    self._all_output_files[sos_variable(x)].append(node)
        self.add_node(node)

    def update_step(self, node, input_targets, depends_targets, output_targets):
        if isinstance(input_targets, sos_targets):
            raise RuntimeError(f"input_targets should not be sos_targets {input_targets}")
        if isinstance(depends_targets, sos_targets):
            raise RuntimeError(f"depends_targets should not be sos_targets {depends_targets}")
        if isinstance(output_targets, sos_targets):
            raise RuntimeError(f"output_targets should not be sos_targets {output_targets}")
        if not isinstance(input_targets, (type(None), Undetermined)):
            for x in input_targets:
                if node not in self._all_dependent_files[x]:
                    self._all_dependent_files[x].append(node)
        if not isinstance(depends_targets, (type(None), Undetermined)):
            for x in depends_targets:
                if node not in self._all_dependent_files[x]:
                    self._all_dependent_files[x].append(node)
        if not isinstance(output_targets, (type(None), Undetermined)):
            for x in output_targets:
                if node not in self._all_output_files[x]:
                    self._all_output_files[x].append(node)

    def find_executable(self):
        '''Find an executable node, which means nodes that has not been completed
        and has no input dependency.'''
        for node in self.nodes():
            # if it has not been executed
            if node._status is None:
                with_dependency = False
                for edge in self.in_edges(node):
                    if edge[0]._status != 'completed':
                        with_dependency = True
                        break
                if not with_dependency:
                    return node
        # if no node could be found, let use try pending ones
        pending_jobs = [x for x in self.nodes() if x._status == 'signature_pending']
        if pending_jobs:
            try:
                notifier = ActivityNotifier(f'Waiting for {len(pending_jobs)} pending job{"s: e.g." if len(pending_jobs) > 1 else ":"} output {short_repr(pending_jobs[0]._signature[0])} with signature file {pending_jobs[0]._signature[1] + "_"}. You can manually remove this lock file if you are certain that no other process is working on the output.')
                while True:
                    for node in pending_jobs:
                        # if it has not been executed
                        lock = fasteners.InterProcessLock(node._signature[1] + '_')
                        if lock.acquire(blocking=False):
                            lock.release()
                            node._status = None
                            return node
                    time.sleep(0.1)
            except Exception as e:
                env.logger.error(e)
            finally:
                notifier.stop()
        return None

    def node_by_id(self, node_uuid):
        for node in self.nodes():
            if node._node_uuid == node_uuid:
                return node
        raise RuntimeError(f'Failed to locate node with UUID {node_uuid}')

    def show_nodes(self):
        for node in self.nodes():
            node.show()
        for edge in self.edges():
            print(edge)

    def circular_dependencies(self):
        try:
            return nx.find_cycle(self)
        except Exception:
            # if there is no cycle, an exception is given
            return []

    def steps_depending_on(self, target, workflow):
        if target in self._all_dependent_files:
            return ' requested by ' + ', '.join(set([workflow.section_by_id(x._step_uuid).step_name() for x in self._all_dependent_files[target]]))
        else:
            return ''

    def pending(self):
        return [x for x in self.nodes() if x._status == 'failed'], [x for x in self.nodes() if x._status is None]

    def dangling(self, targets):
        missing = []
        existing = []
        for x in list(self._all_dependent_files.keys()) + ([] if targets is None else targets):
            if isinstance(x, file_target):
                raise RuntimeError(f'DAG should not contain file_target instance {x}')
            if isinstance(x, sos_targets):
                raise RuntimeError(f'DAG should not contain sos_targets instance {x}')
            if file_target(x).target_exists() if isinstance(x, str) else x.target_exists():
                if x not in self._all_output_files:
                    existing.append(x)
            elif x not in self._all_output_files:
                missing.append(x)
        return missing, existing

    def regenerate_target(self, target):
        if target in self._all_output_files:
            for node in self._all_output_files[target]:
                if node._status == 'completed':
                    if isinstance(target, sos_step):
                        raise RuntimeError(
                            f'Completed target {target} is being re-executed. Please report this bug to SoS developers.')
                    else:
                        env.logger.info(f'Re-running {node._node_id} to generate {target}')
                        node._status = None
            return True
        else:
            # so the signature exists but the step is not in all output files
            # that means the file must have been generated by an auxiliary step
            # but the actual file has been removed.
            # We will have to find this auxiliary step and rerun
            if isinstance(target, str):
                file_target(target).remove_sig()
            else:
                target.remove_sig()
            return False

    def subgraph_from(self, targets):
        '''Trim DAG to keep only nodes that produce targets'''
        # first, find all nodes with targets
        subnodes = []
        for node in self.nodes():
            if not isinstance(node._output_targets, Undetermined) and any(x in node._output_targets for x in targets):
                subnodes.append(node)
        #
        ancestors = set()
        for node in subnodes:
            ancestors |= nx.ancestors(self, node)
        return SoS_DAG(nx.subgraph(self, subnodes + list(ancestors)))

    def build(self, steps):
        '''Connect nodes according to status of targets'''
        # right now we do not worry about status of nodes
        # connecting the output to the input of other nodes
        #
        # NOTE: This is implemented in the least efficient way just for
        # testing. It has to be re-implemented.
        #
        # refer to http://stackoverflow.com/questions/33494376/networkx-add-edges-to-graph-from-node-attributes
        #
        # several cases triggers dependency.
        indexed = [x for x in self.nodes() if x._node_index is not None]
        indexed.sort(key = lambda x: x._node_index)

        for idx, node in enumerate(indexed):
            # 1. if a node changes context (using option alias), all later steps
            # has to rely on it.
            if node._context['__changed_vars__']:
                for later_node in indexed[idx + 1: ]:
                    if node._context['__changed_vars__'] & (later_node._context['__signature_vars__'] | later_node._context['__environ_vars__']):
                        self.add_edge(node, later_node)

            # 2. if the input of a step is undetermined, it has to be executed
            # after all its previous steps.
            if isinstance(node._input_targets, Undetermined) and idx > 0:
                # if there is some input specified, it does not use default
                # input, so the relationship can be further looked before
                if node._input_targets.expr:
                    # if the input is dynamic, has to rely on previous step...
                    if 'dynamic' in node._context['__environ_vars__']:
                        self.add_edge(indexed[idx-1], node)
                    else:
                        # otherwise let us look back.
                        for prev_node in indexed[idx -1 : :-1]:
                            if node._context['__environ_vars__'] & prev_node._context['__changed_vars__']:
                                self.add_edge(prev_node, node)
                else:
                    self.add_edge(indexed[idx-1], node)
        #
        # 3. if the input of a step depends on the output of another step
        for target, in_node in self._all_dependent_files.items():
            for out_node in [y for (x,y) in self._all_output_files.items() if x == target]:
                for i in in_node:
                    for j in out_node:
                        if j != i:
                            self.add_edge(j, i)

    def to_string(self):
        try:
            return nx.drawing.nx_pydot.to_pydot(self).to_string()
        except Exception as e:
            env.logger.warning(f'Failed to call to_pydot: {e}')

    def write_dot(self, filename):
        # write dot, used by tests
        if not filename:
            return
        try:
            nx.drawing.nx_pydot.write_dot(self, filename)
        except Exception as e:
            env.logger.warning(f'Failed to call write_dot: {e}')
