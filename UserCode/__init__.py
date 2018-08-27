import os

usercode_path = __path__[:]


def set_user(username):
    usercode_path.append(username)


def exec_userfile(filename):
    filepath = usercode_path[:]
    filepath.append(filename)
    execfile(os.path.abspath(os.path.join(*filepath)))
