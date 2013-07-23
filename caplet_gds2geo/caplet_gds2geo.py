#!/usr/bin/python

# AUTHOR: Yu-Chung Hsiao
# DATE  : Jan. 31, 2013
# EMAIL : project.caplet@gmail.com
#
# This file is part of CAPLET.
#
# CAPLET is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CAPLET is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with CAPLET.  If not, see <http://www.gnu.org/licenses/>.


import argparse
import sys
import gdsii.library
import gdsii.elements
import collections


def main():
    parser = argparse.ArgumentParser(
        description="Extract geometries from gds2 files")
    parser.add_argument("filename", 
                        help="input gds2 file name", metavar="GDS2_FILE")

    parser.add_argument("-l", "--layer", dest="layerdef_filename",
                        help="layer information", metavar="LAYER_FILE")

    args = parser.parse_args()

    print("LAYER_FILE: " + args.layerdef_filename)
    print("GDSii_FILE: " + args.filename)

    # load the layer definition file
    with open(args.layerdef_filename, 'r') as layerdef:
        (metal, metal_index, via, via_index, via_connect) = parse_layerdef(layerdef)
        
    # set up gds stream
    with open(args.filename, 'rb') as gds_stream:
        lib = gdsii.library.Library.load(gds_stream)


    unit = lib.physical_unit
    struc = lib.pop(0)
    LAYER  = enum('BOTTOM', 'TOP', 'NAME')
    
    # read and sort out structures
    metal_struc = [None]*len(metal_index)
    for i in range(len(metal_struc)):
        metal_struc[i] = collections.deque()

    via_struc = [None]*len(via_index)
    for i in range(len(via_struc)):
        via_struc[i] = collections.deque()

    # only process boundary and path types
    for part in struc:
        if isinstance(part, gdsii.elements.Boundary):
            key = str(part.layer) + '_' + str(part.data_type)
            layer_info = metal.get(key)
            if layer_info != None: # part is in a metal layer
                layer_info = metal.get( layer_info[LAYER.NAME] )
                temp  = layer_info[LAYER.NAME]
                metal_struc[temp].append(part.xy)
                continue

            layer_info = via.get(key)
            if layer_info != None: # part is a via
                layer_info = via.get( layer_info[LAYER.NAME] )
                via_struc[layer_info[LAYER.NAME]].append(part.xy)
                continue
        
        elif isinstance(part, gdsii.elements.Path):
            key = str(part.layer) + '_' + str(part.data_type)
            layer_info = metal.get( metal.get(key)[LAYER.NAME] ) # must be in a metal layer
    
            rects = path2rects(part.xy, part.width)
            for rect in rects:
                metal_struc[layer_info[LAYER.NAME]].append(rect)
           
    # write sorted structures into a geo file
    geo_filename = args.filename.split('.')[0]+'.geo'
    print("OUTPUTFILE: %s" % geo_filename)

    with open(geo_filename, 'w') as output:
        # output metal def
        n_metal = len(metal_index)
        output.write( str( int(n_metal) ) + '\n' )

        for i in range(n_metal):
            output.write( str(i) + ', ' 
                          + str( int(metal[ metal_index[i] ][LAYER.BOTTOM]) ) + ', '
                          + str( int(metal[ metal_index[i] ][LAYER.TOP]) )    + '\n' )

        # output via def
        n_via = len(via_index)
        output.write( str( int(n_via) ) + '\n' )
    
        for i in range(n_via):
            output.write( str(i+n_metal) + ', '
                          + str( int(via[ via_index[i] ][LAYER.BOTTOM]) ) + ', '
                          + str( int(via[ via_index[i] ][LAYER.TOP]) )    + ', '
                          + str( int(metal[ via_connect[ via_index[i] ][0] ][LAYER.NAME]) ) + ', '
                          + str( int(metal[ via_connect[ via_index[i] ][1] ][LAYER.NAME]) ) + '\n' )


        # output metal and struc
        metal_via_struc = [metal_struc, via_struc]
        index = -1
        for each_struc in metal_via_struc:
            for i in range(len(each_struc)):
                index += 1
                output.write( str(index) + '\n' )
                output.write( str(len(each_struc[i])) + '\n' )
                for j in range(len(each_struc[i])):
                    poly = each_struc[i].popleft()
                    output.write( str( len(poly) ) + '\n' )
                    for pnt in poly:
                        output.write( str( int(pnt[0]) ) + ', ' + str( int(pnt[1]) ) + '\n' )

        # complete
        print("Success.")
        

def boundary2point(boundary):
    return [ (boundary[0], boundary[2]),
             (boundary[1], boundary[2]),
             (boundary[1], boundary[3]),
             (boundary[0], boundary[3]),
             (boundary[0], boundary[2])]

def path2rects(xy, width):
    rects = collections.deque()
    width = width/2;

    xy_deque = collections.deque(xy)
    pnt1 = xy_deque.popleft()
    pnt2 = xy_deque[0]
    if pnt1[1]==pnt2[1]: # x-dir
        if pnt1[0]<pnt2[0]: # pnt1 is left to pnt2
            xy_deque.appendleft((pnt1[0]-width, pnt1[1]))
        else: # pnt1 is right to pnt2
            xy_deque.appendleft((pnt1[0]+width, pnt1[1]))
            
    elif pnt1[0]==pnt2[0]: # y-dir
        if pnt1[1]<pnt2[1]: # pnt1 is below pnt2
            xy_deque.appendleft((pnt1[0], pnt1[1]-width))
        else: # pnt1 is above pnt2
            xy_deque.appendleft((pnt1[0], pnt1[1]+width))
            
    pnt1 = xy_deque.pop()
    pnt2 = xy_deque[-1]
    if pnt1[1]==pnt2[1]: # x-dir
        if pnt1[0]<pnt2[0]: # pnt1 is left to pnt2
            xy_deque.append((pnt1[0]+width, pnt1[1]))
        else: # pnt1 is right to pnt2
            xy_deque.append((pnt1[0]-width, pnt1[1]))
            
    elif pnt1[0]==pnt2[0]: # y-dir
        if pnt1[1]<pnt2[1]: # pnt1 is below pnt2
            xy_deque.append((pnt1[0], pnt1[1]+width))
        else: # pnt1 is above pnt2
            xy_deque.append((pnt1[0], pnt1[1]-width))
            
    for i in range(len(xy_deque)-1):
        pnt1 = xy_deque.popleft()
        pnt2 = xy_deque[0]
        
        boundary = [0]*4
        if pnt1[1] == pnt2[1]: # x-dir
            boundary[2] = pnt1[1]-width
            boundary[3] = pnt1[1]+width
            
            if pnt1[0] < pnt2[0]: # pnt1 is left to pnt2
                boundary[0] = pnt1[0]+width
                boundary[1] = pnt2[0]+width
            else: # pnt1 is right to pnt2
                boundary[0] = pnt2[0]-width
                boundary[1] = pnt1[0]-width
                        
        elif pnt1[0] == pnt2[0]: # y-dir
            boundary[0] = pnt1[0]-width
            boundary[1] = pnt1[0]+width
    
            if pnt1[1] < pnt2[1]: # pnt1 is below pnt2
                boundary[2] = pnt1[1]+width
                boundary[3] = pnt2[1]+width
            else: # pnt1 is above pnt2
                boundary[2] = pnt2[1]-width
                boundary[3] = pnt1[1]-width
                     
        rect = boundary2point(boundary)
        rects.append(rect)

    return list(rects)


def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


def parse_layerdef(layerdef):
    STATUS = enum('OUTSIDE', 'METAL', 'VIA')
    LAYER  = enum('BOTTOM', 'TOP')
    METAL  = enum('NAME', 'LAYER', 'DATATYPE', 'HEIGHT', 'THICKNESS')
    VIA    = enum('NAME', 'LAYER', 'DATATYPE', 'BOTTOMNAME', 'TOPNAME')
    # 0: outside blocks
    # 1: in a metal block
    # 2: in a via block

    status = STATUS.OUTSIDE
    # default input unit: um
    unit   = 1e-6
    # defulat output unit: nm
    factor = unit/1e-9

    metal       = dict()
    via         = dict()
    via_connect = dict()

    metal_index = collections.deque()
    via_index = collections.deque()

    lines = layerdef.read()
    blocks = lines.split('}')[0:-1]
    #print('Number of blocks: {0}'.format(str(len(blocks))))

    for block in blocks:
        status = STATUS.OUTSIDE
        lines = block.split('\n')
        for line in lines:
            line = line.split('#')[0] # remove comments
            if line == '': # if empty line
                continue

            if status == STATUS.OUTSIDE:
                if '{' in line:
                    if line.split('{')[0].lower() == 'metal':
                        status = STATUS.METAL
                    elif line.split('{')[0].lower() == 'via':
                        status = STATUS.VIA
                    else:
                        sys.stderr.write('Error: No such a block type: {0}\n'
                                         .format(line.split('{')[0]) )
                        sys.exit(1)

                else:
                    parameter = line.split(':')
                    if parameter[0].lower() == 'unit':
                        try:
                            unit = float(parameter[1])
                        except ValueError:
                            val = parameter[1]
                            if val=='um' or val=='u':
                                unit = 1e-6

                            elif val=='nm' or val=='n':
                                unit = 1e-9
                            elif val=='am' or val=='a':
                                unit = 1e-10
                            else:
                                raise
    
                        factor = unit/1e-9

            elif status == STATUS.METAL:
                try:
                    token = line.split()
                    name = token[METAL.NAME].lower()
                    key  = token[METAL.LAYER]+'_'+token[METAL.DATATYPE]

                    bottom    = int(float(token[METAL.HEIGHT])*factor)
                    thickness = int(float(token[METAL.THICKNESS])*factor)
                    top    = bottom + thickness

                    metal.setdefault(name, [bottom, top, len(metal_index)])
                    metal.setdefault(key,  [bottom, top, name])
                    metal_index.append(name)
                except:
                    print('Format Error: '+line)
                    raise

            elif status == STATUS.VIA:
                try:
                    token = line.split()
                    name = token[VIA.NAME].lower()
                    key  = token[VIA.LAYER]+'_'+token[VIA.DATATYPE]
                    
                    bottom_name = token[VIA.BOTTOMNAME]
                    top_name    = token[VIA.TOPNAME]

                    bottom = metal[bottom_name][LAYER.TOP]
                    top    = metal[top_name][LAYER.BOTTOM]

                    via.setdefault(name, [bottom, top, len(via_index)])
                    via.setdefault(key,  [bottom, top, name])
                    via_connect.setdefault(name, [bottom_name, top_name])
                    via_index.append(name)
                except:
                    print('Format Error: '+line)
                    raise

    return (metal, metal_index, via, via_index, via_connect)

if __name__ == "__main__":
    sys.exit(main())






