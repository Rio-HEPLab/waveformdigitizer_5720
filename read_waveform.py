#!/bin/env python 

import re
import numpy as np
import pandas as pd

if __name__ == '__main__':

    waves = open("wave0.txt", "r") 
    record_length = -1
    event_number = -1
    offset = -1
    arr = None
    p_record_length = re.compile('Record Length: *')
    p_event_number = re.compile('Event Number: *')
    p_offset = re.compile("DC offset \\(DAC\\): *")

    df = pd.DataFrame(columns=['Rec','Evt','DAC','Vals'])

    for idx, line in enumerate(waves): 
        if df.shape[0] >= 1000: break 
	 
        m_record_length = p_record_length.match( line ) 
        m_event_number = p_event_number.match( line ) 
        m_offset = p_offset.match( line ) 
	 
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
        else:     
            if start_array: 
                if arr_entry == -1: 
                    arr = np.zeros( record_length ) 
                    arr_entry = 0 

                val = int( line.rstrip() ) 
                arr[arr_entry] = val 
                arr_entry += 1 
                if arr_entry == (record_length):  
                    print (arr) 
                    df = df.append({'Rec' : record_length, 'Evt' : event_number, 'DAC' : offset, 'Vals': arr}, ignore_index=True) 

    #df.to_csv('output.csv')
    df.to_hdf('output.h5', key='df', mode='w')

