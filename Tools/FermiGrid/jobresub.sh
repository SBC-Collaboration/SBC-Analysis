#################################################################
# condor jobs that were held
#RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499442343"
#for LINE in 19 20 30 34 35 37 43 109
#RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499653396" # 63 runs
#for LINE in 2 3 10 11 12 14 15 17 19 20 21 22 31

#RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499708559" # 24 runs, corrected dytran/acoustic flip
#for LINE in 3 5 6 7 8 12 13 14 15 16 18
#for LINE in 0 1 2 9 10 11 17 19 20 21 22 23 # runs with condor errors

#RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499761912"
#for LINE in 5 6 8 9 10 15 16 17 18 19 21

#RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1499929661" # all runs for Piezo2
#for LINE in 30 32 33 34 35 107 108 109 110 111 120

#RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1500323827"
#for LINE in 4 5 6 7 8 9 10 11 12 14 15 16 17 18 20 21 22 23
#RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1500627181"
#for LINE in 30 32 33 34 37 38 107 108 109 110 111 120 121 194 195 196 197 198 199 200 201 202 203 204 205 206 207 209 210 211 212 213 214 215 216 217 218
#RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/pmtpulse_1503095464"
#for LINE in 2 16 17 18
#RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_SepOct"
#for LINE in 2 3 18 19 20 21 22 23 24 25 26 27 28 59 60 61 62 63 64 66 67 68 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 126 127 128 129 130 131 132 163
RUN_LIST_FILE="/pnfs/coupp/persistent/runlists/SBC-17/list_1510013290"
for LINE in 5 7 13 14 15 16 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 46 47 48 49 54 55 56 57 58 59
do
	LINE=`expr $LINE + 1`
    run=`head -$LINE $RUN_LIST_FILE | tail -1`
    echo $run
    ls -lh /bluearc/storage/SBC-17-data/$run/0/PMTtraces.bin
done
