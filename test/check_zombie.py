import sys
import subprocess
import psutil


def run_process(args):
    proc_before = subprocess.check_output(
        'ps aux | grep -v root', shell=True).decode().splitlines()

    # run the main process
    ret = subprocess.Popen(args)
    ret.wait()

    proc_after = subprocess.check_output(
        'ps aux | grep -v root', shell=True).decode().splitlines()
    pid_before = [
        int(x.split()[1]) for x in proc_before if not 'PID' in x and 'TIME' in x
    ]
    pid_after = [
        int(x.split()[1]) for x in proc_after if not 'PID' in x and 'TIME' in x
    ]
    possible_zombies = [name for pid,name in zip(pid_after, proc_after) if pid not in pid_before and \
        psutil.Process(pid).ppid() not in pid_before]
    if possible_zombies:
        sys.exit('\nNew possible zoombie processes\n{zombies}')
    else:
        print('\nNo possible zombie process is detected')
        sys.exit(ret.returncode)


if __name__ == '__main__':
    if '-h' in sys.argv:
        print('Usage: python check_zombie.py regular command line')
        print(
            '''This command executes the command and lists new processes after the completion of '''
        )
        print(
            '''the command. Processes that are child processes of processes before execution are '''
        )
        print(
            '''excluded. The rest of the processes could be zombie process left by the command, or '''
        )
        print('''new processes created during the execution of the command.''')
    run_process(sys.argv[1:])
