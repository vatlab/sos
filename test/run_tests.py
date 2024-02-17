#!/usr/bin/env python
import argparse
import datetime
import os
import subprocess
import sys

LOGFILE = '.test_results.log'


def get_testcases():
    output = subprocess.check_output(['pytest', '--co'])
    tests = []
    cur_module = ''
    for line in output.decode('utf8').splitlines():
        if line.strip().startswith('<Module'):
            cur_module = line.strip().split()[-1].rstrip('>')
            if not os.path.isfile(cur_module):
                cur_module = f'test/{cur_module}'
            assert os.path.isfile(cur_module)
        if line.strip().startswith('<Function'):
            tests.append(cur_module + '::' + line.strip().split()[-1].rstrip('>'))
    return tests


def run_tests(args, tests, show_output=False):
    failed_tests = []
    if not tests:
        return failed_tests

    def test_failed(test_names, return_code):
        print(f'{" ".join(test_names)}     \x1b[31;1mFAILED\x1b[0m')
        with open(LOGFILE, 'a') as ft:
            ft.write(f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")} {" ".join(test_names)}     FAILED\n')

        if args.exitfirst:
            sys.exit(return_code)
        else:
            failed_tests.extend(test_names)

    try:
        ret = subprocess.run(
            ['pytest'] + list(tests),
            stdout=None if show_output else subprocess.DEVNULL,
            stderr=None if show_output else subprocess.DEVNULL,
            timeout=60 * len(tests))
        if ret.returncode != 0:
            if len(tests) > 1:
                for test in tests:
                    failed_tests.extend(run_tests(args, [test]))
            else:
                test_failed(tests, ret.returncode)
        else:
            with open(LOGFILE, 'a') as log:
                for test in tests:
                    log.write(f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")} {test}     SUCCEED\n')
                    print(f'{test}    \x1b[32;1mPASSED\x1b[0m')
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        if len(tests) > 1:
            for test in tests:
                failed_tests.extend(run_tests(args, [test]))
        else:
            test_failed(tests, 1)
    return failed_tests


if __name__ == '__main__':
    parser = argparse.ArgumentParser('run_tests')
    parser.add_argument('tests', nargs='*', help='Tests that will be executed instead of gathering all tests.')
    parser.add_argument('-b', '--batch', default=5, type=int, help='Group tests')
    parser.add_argument(
        '-l',
        '--lastfailed',
        nargs='?',
        type=int,
        const=0,
        help='''Run only failed tests, default to all. If a number is specified,
        only run the last few failed tests.''')
    parser.add_argument('-x', '--exitfirst', help='Stop when one test fails')
    args = parser.parse_args()

    if args.tests:
        all_tests = args.tests
    else:
        print('Collecting tests')
        all_tests = get_testcases()
        print(f'{len(all_tests)} tests are collected.')

        if args.lastfailed is not None:
            if not os.path.isfile(LOGFILE):
                sys.exit(f'Log file {LOGFILE} does not exists.')
            test_results = {}
            with open(LOGFILE) as fl:
                for line in fl:
                    if not line.strip():
                        continue
                    try:
                        fields = line.split()
                        if len(fields) >= 2:
                            tst = fields[-2]
                            res = fields[-1].strip()
                        else:
                            tst = fields[-1]
                            res = 'FAILED'
                    except Exception:
                        print(f'Invalid log line: {line}')
                    test_results[tst] = res.strip()
            all_tests = [x for x, y in test_results.items() if (y == 'FAILED' and x in all_tests) or x not in test_results]
            # if args.lastfailed != 0:
            #     all_tests = all_tests[-args.lastfailed:]
            print(f'Running {len(all_tests)} failed tests.')

    failed_tests = []
    nbatch = len(all_tests) // args.batch + 1
    for batch in range(nbatch):
        tests = all_tests[batch * args.batch:(batch + 1) * args.batch]
        failed_tests.extend(run_tests(args, tests))

    if failed_tests:
        retried_failed_tests = []
        for test in failed_tests:
            print(f'\n\nRerunning {test}\n')
            retried_failed_tests.extend(run_tests(args, [test], show_output=True))
        #
        failed_tests = retried_failed_tests

    if failed_tests:
        print(f'\n\n{len(failed_tests)} failed tests (logged to {LOGFILE}):\n' + '\n'.join(failed_tests))
    else:
        print(f'All {len(all_tests)} tests complete successfully.')
    sys.exit(0 if not failed_tests else 1)
