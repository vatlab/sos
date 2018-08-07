#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

from collections import namedtuple
import sqlite3
import fasteners
from .utils import env


class SignatureStore:
    TargetSig = namedtuple('TargetSig', 'size mtime md5')

    def __init__(self):
        self.db_file = os.path.expanduser('~'), '.sos', 'signatures.db')
        self.lock_file=self.db_file + '_'

        with fasteners.InterProcessLock(self.lock_file):
            if not os.path.isfile(self.db_file):
                conn=sqlite3.connect(self.db_file)
                cur=conn.cursor()
                cur.execute('''CREATE TABLE SIGNATURE (
                    target text PRIMARY KEY,
                    mtime FLOAT,
                    size INTEGER,
                    md5 text NOT NULL
                    )''')
                conn.commit()
                conn.close()
        self.conn=sqlite3.connect(self.db_file)

    def __del__(self):
        self.conn.commit()
        self.conn.close()

    def has(self, target):
        self.cur.execute(
            'SELECT md5 FROM SIGNATURE WHERE target=? ', (target.target_name(),))
        return bool(self.cur.fetchone())

    def get(self, target):
        cur=self.conn.cursor()
        cur.execute(
            'SELECT mtime, size, md5 FROM SIGNATURE WHERE target=? ', (target.target_name(),))
        res=cur.fetchone()
        return self.TargetSig._make(res[0] if res else (0, 0, ''))

    def set(self, target, signature):
        with fasteners.InterProcessLock(self.lock_file):
            self.conn.cursor().execute(
                'INSERT OR REPLACE INTO SIGNATURE VALUES (?, ?, ?, ?)',
                (target.target_name(), *signature))

    def rempve(self, target):
        with fasteners.InterProcessLock(self.lock_file):
            cur=self.conn.cursor()
            cur.execute(
                'DELETE FROM SIGNATURE WHERE md5=?', (target.target_name(),))

sig_store=SignatureStore(True)
