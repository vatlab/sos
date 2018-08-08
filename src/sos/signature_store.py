#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

from collections import namedtuple
import os
import sqlite3
import fasteners

class TargetSignatures:
    TargetSig = namedtuple('TargetSig', 'mtime size md5')

    def __init__(self):
        self.db_file = os.path.join(os.path.expanduser('~'), '.sos', 'target_signatures.db')
        self.lock_file=self.db_file + '_'
        self._conn = None
        self._pid = None

        with fasteners.InterProcessLock(self.lock_file):
            if not os.path.isfile(self.db_file):
                conn = sqlite3.connect(self.db_file, timeout=20)
                conn.execute('''CREATE TABLE targets (
                    target text PRIMARY KEY,
                    mtime FLOAT,
                    size INTEGER,
                    md5 text NOT NULL
                )''')
                conn.commit()
                conn.close()

    def _get_conn(self):
        # there is a possibility that the _conn is copied with a process
        # and we would better have a fresh conn
        if self._conn is None or self._pid != os.getpid():
            self._conn = sqlite3.connect(self.db_file, timeout=20)
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

    def set(self, target, mtime:float, size:str, md5: str):
        #with fasteners.InterProcessLock(self.lock_file):
        self.conn.cursor().execute(
            'INSERT OR REPLACE INTO targets VALUES (?, ?, ?, ?)',
            (target.target_name(), mtime, size, md5))
        self.conn.commit()

    def remove(self, target):
        #with fasteners.InterProcessLock(self.lock_file):
        cur=self.conn.cursor()
        cur.execute(
                'DELETE FROM targets WHERE md5=?', (target.target_name(),))
        self.conn.commit()

    def clear(self):
        cur=self.conn.cursor()
        cur.execute('DELETE FROM targets')
        cur.execute()
        self.conn.commit()

target_signatures = TargetSignatures()
