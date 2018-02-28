import os

# This is just a helper function I use to run command line stuff
# It is not needed in the workflow that can directly run system commands
def execute_command(command: str, cwd: str=os.getcwd(), shell: bool=True):
    import time
    import sys
    import subprocess
    print("Running process in directory:\n" + cwd)
    print("Command:\n" + command)
    tic = time.time()
    process = subprocess.Popen(command, shell=shell, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True)
    output = []
    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()
        output.append(nextline)
    output = "\n".join(output)
    # output = process.communicate()[0]
    exit_code = process.returncode
    if exit_code == 0:
        return output, sys.stdout, time.time() - tic
    else:
        print("The process ran for " + str(time.time() - tic))
        raise subprocess.CalledProcessError(exit_code, command)