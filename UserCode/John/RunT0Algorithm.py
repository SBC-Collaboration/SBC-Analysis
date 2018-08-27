

from SBCcode.AnalysisModules.EventDealer import ProcessSingleRun2 as psr2











if __name__ == "__main__":
    import os
    event_list_file = "/pnfs/coupp/persistent/runlists/SBC-17/list_1499929661"
    data_dir = "/bluearc/storage/SBC-17-data/"
    with open(event_list_file) as f:
        lines = f.readlines()[-30:]
        events = []
        for line in lines:
            events.append(line.strip())

    for event in events:
        psr2(os.path.join(data_dir, event), dataset="SBC-2017",
             recondir=os.path.join("/pnfs/coupp/persistent/grid_output/SBC-17/output", event),
             process_list = ["acoustic"])
