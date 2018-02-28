import os
import numpy as np

def execute_command(command: str, cwd: str=os.getcwd(), shell: bool=True):
    """
    # This is just a helper function I use to run command line stuff
    # It is not needed in the workflow that can directly run system commands

    :param command: command string
    :param cwd: working directory for the command to run
    :param shell: flag (see subprocess.Popen()
    :return: std output as a string, sys.stdout handle, execution duration in secs
    """
    import time
    import sys
    import subprocess
    print("Running process in directory:\n" + cwd)
    print("Command:\n" + command)
    tic = time.time()
    process = subprocess.Popen(command, shell=shell, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True)
    # TODO: correct it to get std output to output as a string
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



def compute_affine_transform(coords_from, coords_to):
    """
    # Find the affine transformation pomfile -> MRIelectrodes
    :param coords_from: numpy.array with coordinates from source space of size (n_coors, 3)
    :param coords_to:  numpy.array with coordinates from destination space of size (n_coors, 3)
    :return: aff_transform: numpy.array of affine transform (of size 4, 4)
    """
    # Pad the data with ones, so that our transformation can do translations too
    n = coords_from.shape[0]
    pad = lambda x: np.hstack([x, np.ones((x.shape[0], 1))])
    unpad = lambda x: x[:, :-1]
    X = pad(coords_from)
    Y = pad(coords_to)

    # Solve the least squares problem X * A = Y
    # to find our transformation matrix A
    A, res, rank, s = np.linalg.lstsq(X, Y)

    aff_transform = lambda x: unpad(np.dot(pad(x), A))

    print("Target:")
    print(coords_to)
    print("Result:")
    print(aff_transform(coords_from))
    print("Max error:", np.abs(coords_to - aff_transform(coords_from)).max())

    return aff_transform