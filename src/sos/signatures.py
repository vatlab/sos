#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import pickle
import lzma
import sqlite3

from .utils import env

class StepSignatures:
    def __init__(self):
        self.global_db_file = os.path.join(
            os.path.expanduser('~'), '.sos', 'step_signatures.db')
        self.local_db_file = os.path.join(
            env.exec_dir, '.sos', 'step_signatures.db')
        self._global_conn = None
        self._local_conn = None

    def get_conn(self, global_sig=False):
        # there is a possibility that the _conn is copied with a process
        # and we would better have a fresh conn
        if global_sig:
            if self._global_conn is None:
                self._global_conn = sqlite3.connect(
                    self.global_db_file, timeout=60)
                self._global_conn.execute('''CREATE TABLE IF NOT EXISTS steps (
                    step_id text PRIMARY KEY,
                    signature BLOB
                    )''')
            return self._global_conn
        else:
            if self._local_conn is None:
                self._local_conn = sqlite3.connect(
                    self.local_db_file, timeout=60)
                self._local_conn.execute('''CREATE TABLE IF NOT EXISTS steps (
                    step_id text PRIMARY KEY,
                    signature BLOB
                    )''')
            return self._local_conn

    def get(self, step_id: str, global_sig: bool):
        try:
            cur = self.get_conn(global_sig).cursor()
            cur.execute(
                'SELECT signature FROM steps WHERE step_id=? ', (step_id,))
            res = cur.fetchone()
        except sqlite3.DatabaseError as e:
            env.logger.warning(f'Failed to get step signature for step {step_id}: {e}')
            return None
        if res:
            try:
                return pickle.loads(lzma.decompress(res[0]))
            except Exception as e:
                env.logger.warning(
                    f'Failed to load signature for step {step_id}: {e}')
                return None
        else:
            return None

    def set(self, step_id: str, signature: dict, global_sig: bool):
        try:
            conn = self.get_conn(global_sig)
            conn.cursor().execute(
                'INSERT OR REPLACE INTO steps VALUES (?, ?)',
                (step_id, lzma.compress(pickle.dumps(signature))))
            conn.commit()
        except sqlite3.DatabaseError as e:
            env.logger.warning(f'Failed to set step signature for step {step_id}: {e}')

    def _num_records(self, cur):
        try:
            cur.execute('SELECT COUNT(rowid) FROM steps')
            return cur.fetchone()[0]
        except sqlite3.DatabaseError as e:
            env.logger.warning(f'Failed to number of records for steps: {e}')
            return 0

    def remove_many(self, steps: list, global_sig: bool):
        try:
            conn = self.get_conn(global_sig)
            cur = conn.cursor()
            cnt = self._num_records(cur)
            cur.executemany(
                    'DELETE FROM steps WHERE step_id=?',
                    [(x,) for x in steps])
            conn.commit()
            return cnt - self._num_records(cur)
        except sqlite3.DatabaseError as e:
            env.logger.warning(f'Failed to remove signature for {len(steps)} substeps: {e}')
            return 0

    def clear(self, global_sig: bool):
        try:
            conn = self.get_conn(global_sig)
            conn.execute('DELETE FROM steps')
            conn.commit()
        except sqlite3.DatabaseError as e:
            env.logger.warning(f'Failed to clear step signature database: {e}')

    def close(self):
        self.get_conn(True).close()
        self.get_conn(False).close()

class WorkflowSignatures(object):
    def __init__(self):
        self.db_file = os.path.join(
            env.exec_dir, '.sos', 'workflow_signatures.db')
        self._conn = None

    def _get_conn(self):
        # there is a possibility that the _conn is copied with a process
        # and we would better have a fresh conn
        if self._conn is None:
            self._conn = sqlite3.connect(self.db_file, timeout=60)
            self._conn.execute('''CREATE TABLE IF NOT EXISTS workflows (
                    master_id text,
                    entry_type text,
                    id text,
                    item text
            )''')
        return self._conn

    conn = property(_get_conn)

    def write(self, entry_type: str, id: str, item: str):
        try:
            self.conn.execute('INSERT INTO workflows VALUES (?, ?, ?, ?)',
                      (env.config["master_id"], entry_type, id, item))
            self.conn.commit()
        except sqlite3.DatabaseError as e:
            env.logger.warning(f'Failed to write workflow signature of type {entry_type} and id {id}: {e}')
            return None

    def records(self, workflow_id):
        try:
            cur = self.conn.cursor()
            cur.execute(
                'SELECT entry_type, id, item FROM workflows WHERE master_id = ?', (workflow_id,))
            return cur.fetchall() if cur else []
        except sqlite3.DatabaseError as e:
            env.logger.warning(f'Failed to get records of workflow {workflow_id}: {e}')
            return []

    def workflows(self):
        try:
            cur = self.conn.cursor()
            cur.execute('SELECT DISTINCT master_id FROM workflows')
            return [x[0] for x in cur.fetchall()]
        except sqlite3.DatabaseError as e:
            env.logger.warning(f'Failed to get workflows from signature database: {e}')
            return []

    def tasks(self):
        try:
            cur = self.conn.cursor()
            cur.execute('SELECT DISTINCT id FROM workflows WHERE entry_type = "task"')
            return [x[0] for x in cur.fetchall()]
        except sqlite3.DatabaseError as e:
            env.logger.warning(f'Failed to get tasks from signature database: {e}')
            return []

    def files(self):
        '''Listing files related to workflows related to current directory'''
        try:
            cur = self.conn.cursor()
            cur.execute('SELECT id, item FROM workflows WHERE entry_type = "tracked_files"')
            return [(x[0], eval(x[1])) for x in cur.fetchall()]
        except sqlite3.DatabaseError as e:
            env.logger.warning(f'Failed to get files from signature database: {e}')
            return []

    def placeholders(self, workflow_id = None):
        try:
            cur = self.conn.cursor()
            if workflow_id is None:
                cur.execute('SELECT item FROM workflows WHERE entry_type = "placeholder"')
            else:
                cur.execute(f'SELECT item FROM workflows WHERE entry_type = "placeholder" AND master_id = "{workflow_id}"')
            return [x[0] for x in cur.fetchall()]
        except sqlite3.DatabaseError as e:
            env.logger.warning(f'Failed to get placeholders from signature database: {e}')
            return []

    def clear(self):
        try:
            self.conn.execute(
                f'DELETE FROM workflows WHERE master_id = ?', (env.config["master_id"],))
            self.conn.commit()
        except sqlite3.DatabaseError as e:
            env.logger.warning(f'Failed to clear workflow database: {e}')
            return []

    def close(self):
        self.conn.close()
