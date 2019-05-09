#!/usr/bin/env python
import sys
sys.argv.append( '-b-' )
import ROOT as r
import json
r.gROOT.SetBatch(True)
sys.argv.remove( '-b-' )

from ROOT import TFile
from pprint import pprint

def getTripletFromFile(fName):
    f = TFile.Open(fName)
    # f.Print()
    t = f.Get('limit')
    # t.Scan('*')
    vals = []
    for e in t:
        if e.quantileExpected == -1:
            continue
        vals.append( (e.quantileExpected, e.limit) )

    # print vals
    return UncertaintiesFromTriplet(vals)

def UncertaintiesFromTriplet(triplet):
    #This assumes that the values are ordered as central(0),lower(1),upper(2)
    #This could be found from the quantileExpected value, but I was lazy
    up = triplet[2][1] - triplet[0][1]
    dn = triplet[1][1] - triplet[0][1]
    # print dn, up
    return (dn, up)

from math import sqrt
def quadDiff(unc1, unc2):
    dn = sqrt( abs(unc1[0]**2 - unc2[0]**2 ))
    up = sqrt( abs(unc1[1]**2 - unc2[1]**2) )
    return (-dn, up)

uncs = []
for f in (sys.argv[1],sys.argv[2]):
    uncs.append( getTripletFromFile(f) )

diff = quadDiff(uncs[0],uncs[1])

#print diff
fileIn = sys.argv[4]

json_file = open(fileIn, 'r')
data = json.load(json_file)
nom_afb = data['POIs'][0]['fit'][1]
my_param = {}
my_param['Afb'] = [nom_afb + diff[0], nom_afb, nom_afb + diff[1]]
my_param['fit'] = [1.0,1.0,1.0]
my_param['groups'] = []
my_param['impact_Afb'] = max(abs(diff[0]), abs(diff[1]))
my_param['name'] = sys.argv[3]
my_param['prefit'] = [1.0,1.0,1.0]
my_param["type"]= "Unconstrained"
data['params'].append(my_param)
json_file.close()
out_file = open(fileIn, 'w')
json.dump(data, out_file, separators=(', \n', ': \n'))


