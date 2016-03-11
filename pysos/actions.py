
import subprocess
import tempfile
import pipes

from .utils import env


#
# A decoration function that allows SoS to replace all SoS actions
# with a null action.
#
def SoS_Action(func):
    def action(*args, **kwargs):
        if env.run_mode == 'dryrun':
            return 0
        else:
            return func(*args, **kwargs)
    return action


class SoS_ExecuteScript:
    def __init__(self, script, interpreter, suffix):
        self.script = script
        self.interpreter = interpreter
        self.script_file = tempfile.NamedTemporaryFile(mode='w+t', suffix=suffix, delete=False).name
        with open(self.script_file, 'w') as script_file:
            script_file.write(self.script)

    def run(self):
        if '{}' in self.interpreter:
            cmd = self.interpreter.replace('{}', pipes.quote(self.script_file))
        else:
            cmd = self.interpreter + ' ' + pipes.quote(self.script_file) 
        env.logger.info('Running ``{}``'.format(cmd))
        ret = subprocess.call(cmd, shell=True)
        if ret != 0:
            raise RuntimeError('Failed to execute script')


class SoS_ExecuteRScript(SoS_ExecuteScript):
    '''SoS_Execute in-line R script using Rscript as interpreter. Please
    check action SoS_ExecuteScript for more details.
    '''
    def __init__(self, script=''):
        SoS_ExecuteScript.__init__(self, script=script, interpreter='Rscript',
            suffix='.R')

class SoS_ExecuteShellScript(SoS_ExecuteScript):
    '''SoS_Execute in-line shell script using bash as interpreter. Please
    check action SoS_ExecuteScript for more details.
    '''
    def __init__(self, script=''):
        SoS_ExecuteScript.__init__(self, script=script, interpreter='bash', 
            suffix='.sh')

class SoS_ExecutePythonScript(SoS_ExecuteScript):
    '''SoS_Execute in-line python script using python as interpreter. Please
    check action SoS_ExecuteScript for more details.
    '''
    def __init__(self, script=''):
        SoS_ExecuteScript.__init__(self, script=script, interpreter='python',
            suffix='.py')


class SoS_ExecutePython3Script(SoS_ExecuteScript):
    '''SoS_Execute in-line python script using python as interpreter. Please
    check action SoS_ExecuteScript for more details.
    '''
    def __init__(self, script=''):
        SoS_ExecuteScript.__init__(self, script=script, interpreter='python3',
            suffix='.py')

@SoS_Action
def run(script):
    SoS_ExecuteShellScript(script).run()

@SoS_Action
def python(script):
    SoS_ExecutePythonScript(script).run()

@SoS_Action
def python3(script):
    SoS_ExecutePython3Script(script).run()

@SoS_Action
def R(script):
    SoS_ExecuteRScript(script).run()

