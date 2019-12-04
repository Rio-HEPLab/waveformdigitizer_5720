#!/bin/env python 

import sys
import re
import numpy as np
import pandas as pd

def main():
    
    filename = "doubleChannel.txt"
    maxEvents = 10000

    #for arg in sys.argv[1:]:
    #    if isinstance(arg, str):
    #        print("Reading File " + arg)
    #        filename = arg
    #    if isinstance(arg, int):
    #        print("Number of events " + arg)
    #        maxEvents = arg

    args = sys.argv[1:]
    if len(args) >= 1:
        filename = args[0]     
    if len(args) >= 2:
        maxEvents = int(args[1])

    print("Number of events {:d}".format(maxEvents))
    print("Reading File " + filename)

    waves = open(filename, "r") 
    record_length = -1
    event_number = -1
    offset = -1
    arr = []
    channel = -1
    p_record_length = re.compile('Record Length: *')
    p_event_number = re.compile('Event Number: *')
    p_offset = re.compile("DC offset \\(DAC\\): *")
    p_channel = re.compile('Channel:*')

    df = pd.DataFrame(columns=['Rec','Evt','DAC','Vals'])

    for idx, line in enumerate(waves): 
        if maxEvents >= 0 and df.shape[0] >= maxEvents: break 
	 
        m_record_length = p_record_length.match( line ) 
        m_event_number = p_event_number.match( line ) 
        m_offset = p_offset.match( line )
        m_channel = p_channel.match( line )
	 
        if m_record_length: 
            record_length = int( line[m_record_length.end():] )  
            print (line, record_length) 
            start_array = False 
            arr_entry = -1 
        elif m_event_number:     
            event_number = int( line[m_event_number.end():] ) 
            print (line, event_number) 
        elif m_offset:      
            offset = int( line[m_offset.end():], 0 ) 
            print (line, offset) 
            start_array = True 
        elif m_channel:
            channel = int(line[m_channel.end():])
            print(line, channel)
        else:     
            if start_array: 
                if arr_entry == -1: 
                    arr.append(np.zeros( record_length ))
                    arr_entry = 0 

                val = int( line.rstrip() ) 
                arr[channel][arr_entry] = val 
                arr_entry += 1 
                if arr_entry == (record_length) and channel > 0:  
                    print (arr)
                    df = df.append({'Rec' : record_length, 'Evt' : event_number, 'DAC' : offset, 'Vals': arr}, ignore_index=True) 
                    arr = None

    #df.to_csv('output.csv')
    #df.to_hdf('output.h5', key='df', mode='w')
    #df.to_hdf('output-wave1-SingleTrig-V2-40k.h5', key='df', mode='w')
    file_str = ( filename.split('/')[-1] ).split('.')[0]
    df.to_hdf('output-' + file_str + '.h5', key='df', mode='w')

if __name__== '__main__':
    main()
