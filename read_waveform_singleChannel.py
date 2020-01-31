#!/bin/env python

import argparse
import re
import numpy as np
import h5py

def main():

    #maxEvents = 10000
    #for arg in sys.argv[1:]:
    #    if isinstance(arg, str):
    #        print("Reading File " + arg)
    #        filename = arg
    #    if isinstance(arg, int):
    #        print("Number of events " + arg)
    #        maxEvents = arg

    parser = argparse.ArgumentParser(description = 'Converte saida do Wavedump para hdf5')
    parser.add_argument('fileName', help = 'Arquivo de entrada' )
    parser.add_argument('-n', '--events', dest = 'events', type = int, required = False, default = 10000, help = 'Numero de eventos' )
    #parser.add_argument('-f', action = 'store', dest = 'txt', required = True, help = 'Saida do Wavedump' )
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'debug', required = False, help = 'Flag de debug' )
    args = parser.parse_args()

    fileName = args.fileName
    events = args.events
    print( "Number of events {:d}".format( events ) )
    print( "Reading File " + fileName )

    file_str = ( fileName.split('/')[-1] ).split('.')[0]
    with h5py.File('output-' + file_str + '.h5', 'w') as f:

        waves = open(fileName, "r")
        record_length = -1
        event_number = -1
        offset = -1
        arr = []
        p_record_length = re.compile('Record Length: *')
        p_event_number = re.compile('Event Number: *')
        p_offset = re.compile("DC offset \\(DAC\\): *")

        data = []

        i = 0
        for idx, line in enumerate(waves):
            if events >= 0 and i >= events: break

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
                        i += 1
                        print('Iteration', i)
                        data.append(arr)
                        arr = None

        dset = f.create_dataset('Vals', data=data)

if __name__== '__main__':
    main()
