import os
import subprocess


def pathjoin(*parts):
    """Because os.pathjoin does unintuitive stuff like
    os.path.join("/home/foo", "/bar/") returning '/bar/'
    """
    return "/".join(map(lambda part: str(part).strip("/"), parts))


def run_cmd(cmd, throw_on_error=True, env=None, stream_output=False, **kwargs):
    """Runs a command as a child process.

    A convenience wrapper for running a command from a Python script.
    Keyword arguments:
    cmd -- the command to run, as a list of strings
    throw_on_error -- if true, raises an Exception if the exit code of the program is nonzero
    env -- additional environment variables to be defined when running the child process
    stream_output -- if true, does not capture standard output and error; if false, captures these
      streams and returns them

    Note on the return value: If stream_output is true, then only the exit code is returned. If
    stream_output is false, then a tuple of the exit code, standard output and standard error is
    returned.
    """
    cmd_env = os.environ.copy()
    if env:
        cmd_env.update(env)

    if stream_output:
        child = subprocess.Popen(cmd, env=cmd_env, **kwargs)
        exit_code = child.wait()
        if throw_on_error and exit_code is not 0:
            raise Exception("Non-zero exitcode: %s" % (exit_code))
        return exit_code
    else:
        child = subprocess.Popen(
            cmd, env=cmd_env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
        (stdout, stderr) = child.communicate()
        exit_code = child.wait()
        if throw_on_error and exit_code is not 0:
            raise Exception("Non-zero exitcode: %s\n\nSTDOUT:\n%s\n\nSTDERR:%s" % (exit_code,
                                                                                   stdout, stderr))
        return (exit_code, stdout, stderr)


def run_streaming_cmd(cmd):
    """Run a command and return a generator that yields each line.

    If the process completes with exit code 0, the generator will complete successfully, otherwise
    we thrown an excption containing the error code.
    """
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        retcode = process.poll()  # returns None while subprocess is running
        line = process.stdout.readline()
        yield line
        if retcode is not None:
            if retcode == 0:
                break
            else:
                raise Exception("Non-zero ret-code: %s" % retcode)


def run_check_output(cmd, shell=False):
    """Run command and return its output as a byte string (NOT return code)."""
    # https://docs.python.org/2/library/subprocess.html#frequently-used-arguments
    if shell:
        if isinstance(cmd, list):
            command = " ".join(cmd)
        else:
            command = cmd
    else:
        if isinstance(cmd, list):
            command = cmd
        else:
            if '"' in cmd or "'" in cmd:
                raise NotImplementedError(
                    "Parsing of single string input containing quotes is not supported. "
                    "Pass command in as list of strings instead.")
            command = cmd.split(" ")
    try:
        output = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=shell)
        if type(output) is bytes:
            return output.decode("utf-8")
        else:
            return output
    except subprocess.CalledProcessError as e:
        print("Command {} failed".format(command))
        print("Error message: {}".format(str(e)))
        print("Return code: {}".format(e.returncode))
        print("Output: {}".format(e.output))
        raise e


def which(program):
    """Returns the path to a program on the caller's system, or None if the program is not present.
    Can be replaced by shutil.which() in python 3.3.
    From http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python"""

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    # If program is a relative or absolute path, check the file at that location.
    # Otherwise, append it to each directory in the PATH environment variable and check whether an
    # executable file exists at that path..
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None