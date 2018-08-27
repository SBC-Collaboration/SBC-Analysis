import numpy as np


def EventAnalysis(ev):
    default_output = dict(Event_run_type=np.int32(-1),
                          Event_trigger_main=np.int32(-1),
                          Event_trigger_cameras=np.int32(-1),
                          Event_trigger_PLC=np.int32(-1),
                          Event_trigger_slowDAQ=np.int32(-1),
                          Event_timestamp=np.float64(-1),
                          Event_mstick=np.int64(-1),
                          Event_Pset=np.float64(-1),
                          Event_livetime=np.float64(-1)
                          )
    try:
        if not ev['event']['loaded']:
            return default_output
        out = default_output
        for k in ev['event']:
            if not (k == 'loaded'):
                out['Event_' + k] = ev['event'][k]
        return out
    except:
        return default_output
