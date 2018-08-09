#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

from collections import namedtuple
import os
import pickle
import lzma
import sqlite3

from .utils import env


class TargetSignatures:
    TargetSig = namedtuple('TargetSig', 'mtime size md5')

    def __init__(self):
        self.db_file = os.path.join(os.path.expanduser(
            '~'), '.sos', 'target_signatures.db')
        self._conn = None
        self._pid = None

    def _get_conn(self):
        # there is a possibility that the _conn is copied with a process
        # and we would better have a fresh conn
        if self._conn is None or self._pid != os.getpid():
            self._conn = sqlite3.connect(self.db_file, timeout=20)
            self._conn.execute('''CREATE TABLE IF NOT EXISTS targets (
                target text PRIMARY KEY,
                mtime FLOAT,
                size INTEGER,
                md5 text NOT NULL
                )''')
            self._pid = os.getpid()
        return self._conn

    conn = property(_get_conn)

    def _list_all(self):
        cur = self.conn.cursor()
        cur.execute('SELECT * FROM targets;')
        for rec in cur.fetchall():
            print(self.TargetSig._make(rec))

    def get(self, target):
        cur = self.conn.cursor()
        cur.execute(
            'SELECT mtime, size, md5 FROM targets WHERE target=? ', (target.target_name(),))
        res = cur.fetchone()
        return self.TargetSig._make(res) if res else None

    def set(self, target, mtime: float, size: str, md5: str):
        self.conn.cursor().execute(
            'INSERT OR REPLACE INTO targets VALUES (?, ?, ?, ?)',
            (target.target_name(), mtime, size, md5))
        self.conn.commit()

    def remove(self, target):
        cur = self.conn.cursor()
        cur.execute(
            'DELETE FROM targets WHERE target=?', (target.target_name(),))
        self.conn.commit()

    def clear(self):
        cur = self.conn.cursor()
        cur.execute('DELETE FROM targets')
        self.conn.commit()


class StepSignatures:
    def __init__(self):
        self.global_db_file = os.path.join(
            os.path.expanduser('~'), '.sos', 'step_signatures.db')
        self.local_db_file = os.path.join(
            env.exec_dir, '.sos', 'step_signatures.db')
        self._global_conn = None
        self._local_conn = None
        self._pid = None

    def get_conn(self, global_sig=False):
        # there is a possibility that the _conn is copied with a process
        # and we would better have a fresh conn
        if global_sig:
            if self._global_conn is None or self._pid != os.getpid():
                self._global_conn = sqlite3.connect(
                    self.global_db_file, timeout=20)
                self._global_conn.execute('''CREATE TABLE IF NOT EXISTS steps (
                    step_id text PRIMARY KEY,
                    signature BLOB
                    )''')
                self._pid = os.getpid()
            return self._global_conn
        else:
            if self._local_conn is None or self._pid != os.getpid():
                self._local_conn = sqlite3.connect(
                    self.local_db_file, timeout=20)
                self._local_conn.execute('''CREATE TABLE IF NOT EXISTS steps (
                    step_id text PRIMARY KEY,
                    signature BLOB
                    )''')
                self._pid = os.getpid()
            return self._local_conn

    def get(self, step_id: str, global_sig: bool):
        cur = self.get_conn(global_sig).cursor()
        cur.execute(
            'SELECT signature FROM steps WHERE step_id=? ', (step_id,))
        res = cur.fetchone()
        if res:
            try:
                return pickle.loads(lzma.decompress(res[0]))
            except Exception as e:
                env.logger.warning(
                    f'Failed to load signature for step {step_id}')
                return None
        else:
            return None

    def set(self, step_id: str, signature: dict, global_sig: bool):
        # with fasteners.InterProcessLock(self.lock_file):
        conn = self.get_conn(global_sig)
        conn.cursor().execute(
            'INSERT OR REPLACE INTO steps VALUES (?, ?)',
            (step_id, lzma.compress(pickle.dumps(signature))))
        conn.commit()

    # def remove(self, step_id: str, global_sig:bool):
    #     #with fasteners.InterProcessLock(self.lock_file):
    #     conn = self.get_conn(global_sig)
    #     cur=self.conn.cursor()
    #     cur.execute(
    #             'DELETE FROM steps WHERE step_id=?', (step_id,))
    #     conn.commit()
    #
    def clear(self, global_sig: bool):
        conn = self.get_conn(global_sig)
        conn.execute('DELETE FROM steps')
        conn.commit()


class WorkflowSignatures(object):
    def __init__(self):
        self.db_file = os.path.join(
            env.exec_dir, '.sos', 'workflow_signatures.db')
        self._conn = None
        self._pid = None

    def _get_conn(self):
        # there is a possibility that the _conn is copied with a process
        # and we would better have a fresh conn
        if self._conn is None or self._pid != os.getpid():
            self._conn = sqlite3.connect(self.db_file, timeout=20)
            self._conn.execute('''CREATE TABLE IF NOT EXISTS workflows (
                    master_id text,
                    entry_type text,
                    id text,
                    item text
            )''')
            self._pid = os.getpid()
        return self._conn

    conn = property(_get_conn)

    def write(self, entry_type: str, id: str, item: str):
        self.conn.execute('INSERT INTO workflows VALUES (?, ?, ?, ?)',
                          (env.config["master_id"], entry_type, id, item))
        self.conn.commit()

    def records(self, workflow_id):
        cur = self.conn.cursor()
        cur.execute(
            'SELECT entry_type, id, item FROM workflows WHERE master_id = ?', (workflow_id,))
        return cur if cur else []

    def workflows(self):
        cur = self.conn.cursor()
        cur.execute('SELECT DISTINCT master_id FROM workflows')
        return [x[0] for x in cur.fetchall()]

    def clear(self):
        self.conn.execute(
            f'DELETE FROM workflows WHERE master_id = ?', (env.config["master_id"],))
        self.conn.commit()


target_signatures = TargetSignatures()
step_signatures = StepSignatures()
workflow_signatures = WorkflowSignatures()
