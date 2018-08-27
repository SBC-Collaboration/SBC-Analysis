from .DataHandling.ReadBinary import ReadBlock as read_bin
from .DataHandling.ReadText import ReadFile as read_txt
from .DataHandling.GetSBCEvent import GetEvent as get_event
from .AnalysisModules.EventDealer import ProcessSingleRun2 as psr2
__all__ = ["DataHandling"]


def add_usercode_to_path(username):
    import sys
    import os.path
    userpath = __path__[:]
    userpath.append('UserCode')
    userpath.append(username)
    userpath = os.path.abspath(os.path.join(*userpath))
    os.path.walk(userpath, lambda a, b, c: sys.path.append(b), 0)
