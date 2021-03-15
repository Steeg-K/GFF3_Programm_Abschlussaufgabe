# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 20:12:30 2021

@author: morit
"""
#test

import logging
import os.path
import matplotlib.pyplot as plt
import argparse
# Ursprungsdatei für die Weiterverarbeitung vorbereiten
logging.basicConfig(level=20)
def speichern(file):
    
    # Zeilen der gff3 in Liste speichern
    liste = []
    buchstabe = True

    while buchstabe:
        for line in file:
            line = line.strip()
            if not line.startswith(">"):
                indicec = line.split('\t')
                liste.append(indicec)
            else:
                buchstabe = False
    
    # Für FASTA-Aufbereitung
    file.seek(0)
    plasmide = {}
    anzahl = []
    zähler = 0
    seq = ""
    ohnev = ""

    for line in file:
        
        #guckt, ob eine neue separate Sequenz anfängt
        if line.endswith("\n") and line.startswith(">"):
            zeile = line.strip()
            
            #String mit Sequenz wird zusammen mit Sequenznamen ins dict "plasmide" eingetragen
            plasmide[ohnev] = seq.lower()
            
            #String reset für nächste Sequenz
            seq = ""
            
            
            ohnev = zeile.split(">")[1]
            zähler += 1 
            anzahl.append(line)
        
        #speichert alle Zeilen, die nicht der Sequenzname sind, ab
        if line.endswith("\n") and zähler > 0:
            if line.startswith(">"):
                pass
            else:
                zeile  = line.strip()
                
                # Zähler zählt mit, wie viele separate Sequenzen da sind. Wenn eine neue beginnt, wird der string seq resettet, um die nächste Sequenz abzuspeichern
                # seq wird nur erweitert, bis eine neue seq beginnt (der zähler ist dann nicht mehr gleich anzhal, weil dieses um einen Eintrag erweitert wurde)
                if zähler == len(anzahl):
                     
                    #speichert jede Zeile aneinandergehängt in einen String
                    seq += zeile
    
    #weil am Ende des files zwar ohnev und seq mit den Daten der letzten Sequenz gefüllt sind, aber der loop vor dem letzten Eintrag ins dict aufhört, kommt der Eintrag hier
    plasmide[ohnev] = seq.lower()
    
    #entfernen des leeren Eintrags am Anfang, der durch die Programmstruktur bedingt war
    plasmide.pop("")  
    
    
    return liste,plasmide


#suchen nach Merkmalen und ihrer Anzahl
def suchen(file):
    zähler = {}
    for line in file:
        # alle Zeilen überspringe, die nicht im Tabellenformat der gff3-Datei sind mit 9 verschiedenen Fields
        if len(line) == 9:
            
            #Anlegen eines dict mit allen Merkmalen und deren Häufigkeit in der Datei
            try:
                zähler[line[2]] += 1
            except KeyError:
                zähler[line[2]] = 1
            
            #letzte Spalte mit den Attributen wird gesplittet, sodass die Attribute einzeln in einer Liste stehen
            attributes = line[8].split(";")
            
            #Zählen der CDS mit hypothetical protein als Namen
            if "Name=hypothetical protein" in attributes and line[2] == "CDS":
                try:
                    zähler["CDS mit hypothetical protein"] +=1
                except KeyError:
                    zähler["CDS mit hypothetical protein"] = 1
            
            #Vorbereitung zum Filtern nach dem Vorkommen einzelner Attributsarten, ohne den Wert dazu zu kennen
            aufteilung = ""
            for index,inhalt in enumerate(attributes):
                aufteilung += "\n " + str(attributes[index])
            
            #Zählen der Merkmale mit annotierem "gene"
            if "gene" in aufteilung:
                try:
                    zähler["Merkmale mit annotiertem gene"] +=1
                except KeyError:
                    zähler["Merkmale mit annotiertem gene"] = 1
            
            #Zählen der MErkmale mit COG category
            if "Dbxref" in aufteilung:
                    try:
                        aufteilung.find("COG:")
                        zähler["Merkmale mit COG category"] +=1
                    except KeyError:
                        zähler["Merkmale mit COG category"] = 1
                    except :
                        pass
                
    print("\nannotierte Sequenzabschnitte samt Anzahl:")          
    print(zähler)              
    return zähler


def übersicht(file, filename):
    ohneendung = filename.split(".")[0]
    outfile = open("Übersicht_"+ ohneendung + ".csv","w")
    zusammenfassung = "#Name,Feature type,Start,Stop,Gene name,Product,Strand\n"
    
    for line in file:      
        if len(line) == 9:
           attributes = line[8].split(";")
           name = ""
           genname = ""
           produkt = ""


            #Suchen der Attribute, die für die Übersicht relevant sind
           for index,inhalt in enumerate(attributes):
                aufteilung =  str(attributes[index])
                
                if aufteilung.startswith("Name="):
                    #Umgehen der Problematik, die durch Kommata in den Namen ensteht (damit nicht zu viele Zellen erzeugt werden dann)
                    if aufteilung.find(",") < 0:
                        ohnegleich = aufteilung.split("=")[1]
                        name = ohnegleich
                    else: 
                        ohnegleich = aufteilung.split("=")[1]
                        ohnegleich = "\"" + ohnegleich + "\""
                        name = ohnegleich
                if aufteilung.startswith("product="):
                     if aufteilung.find(",") < 0:
                        ohnegleich = aufteilung.split("=")[1]
                        produkt = ohnegleich
                     else: 
                        ohnegleich = aufteilung.split("=")[1]
                        ohnegleich = "\"" + ohnegleich + "\""
                        produkt = ohnegleich
                if aufteilung.startswith("gene="):
                    ohnegleich = aufteilung.split("=")[1]
                    genname = ohnegleich
           
          #Erzeugen der Zeile mit den ganzen Basisinformationen eines MErkmals
           zusammenfassung += name + "," 
           zusammenfassung += line[2] + "," 
           zusammenfassung += line[3] + "," 
           zusammenfassung += line[4] + "," 
           zusammenfassung += genname + "," 
           zusammenfassung += produkt + "," 
           zusammenfassung += line[6] + "\n" 

    
    #Erzeugen der Datei
    outfile.writelines(zusammenfassung)
    outfile.close()
    
    logging.info("Übersicht als CSV erstellt\n")
            
def codoncount(seq, gff3):
    
    anzahl = {}
    anteil = {}
    insgesamt = 0
    startinsg = 0
    startcodons = {}
    startanteil = {}
    
    for key,value in seq.items():
        
       for element in gff3:
           if len(element) == 9 and key == element[0] and element[2] == "CDS": 
               start = int(element[3])
               stop = int(element[4])
               
               if stop > len(value):
                   print("\n")
                   logging.warning("\nannotierte CDS "+ key +"_überschreitet Sequenzlänge und wird deswegen nicht ausgewertet")
               else:
                       
                   if element[6] == "+":
                           
                       
                       for char in range(start-1,stop-3,3):
                            
                               if char == start-1:
                                   starttriplet = value[char] + value[char+1] + value[char+2]
                                   startcodon = starttriplet.replace("t","u")
                                   try:
                                       startcodons[startcodon] += 1
                                   except KeyError:
                                       startcodons[startcodon] = 1
                                       
                               triplet = value[char] + value[char+1] + value[char+2]
                               
                               codon = triplet.replace("t","u")
                               
                               try:
                                   anzahl[codon] += 1
                               except KeyError:
                                   anzahl[codon] = 1
                  
                   if element[6]=="-":
                        komp_seq = komplementieren(value,start,stop)
                       
                        for char in range(0,len(komp_seq)-2,3):
                            
                               if char == 0:
                                   starttriplet = komp_seq[char] + komp_seq[char+1] + komp_seq[char+2]
                                   startcodon = starttriplet.replace("t","u")
                                   try:
                                       startcodons[startcodon] += 1
                                   except KeyError:
                                       startcodons[startcodon] = 1
                                       
                               triplet = komp_seq[char] + komp_seq[char+1] + komp_seq[char+2]
                               
                               codon = triplet.replace("t","u")
                               
                               try:
                                   anzahl[codon] += 1
                               except KeyError:
                                   anzahl[codon] = 1
                   
    for key, value in anzahl.items():
        insgesamt += value
        
    for key, value in startcodons.items():
        startinsg += value
    
    for key, value in anzahl.items():
        anteil[key] = round(value/insgesamt*100,2) 
        
    for key,value in startcodons.items():
        startanteil[key] = round(value/startinsg*100,2)
    
    print("\nCodonusage:")
    print(anzahl)
    print("\nProzentualer Anteil der genutzten Codons:")
    print(anteil)
    print("\nAnzahl der genutzten Startcodons:")
    print(startcodons)
    print("\nAnteil der genutzten Startcodons von der Gesamtheit:") 
    print(startanteil)
    
    return anzahl,startcodons   

def genextraktion(file):
    
    plasmide = []
    merkmalsammlung = []
    dictio = {}
    
    
    for line in file:
        
        #wenn das format erkannt wird, in dem cds annotiert wurden
        if len(line) == 9:
            if line[2] == "CDS":
                
                
                #wenn ein neues plasmid/seq gelesen wird
                if line[0] not in plasmide:
                    
                    #letztes plasmid im dict abspeicerhn mit informationen der einzelnen cds aus merkmalsammlung
                    try:
                        dictio[merkmalsammlung[0]] = merkmalsammlung[1:]
                        merkmalsammlung.clear()
                    except:
                        pass
                    
                    
                    #neues plasmid merken, anfangen mit dem abspeicehrn der zeileninformationen für erste cds
                    plasmide.append(line[0])
                    merkmalsammlung.append(line[0])
                    start = line[3]
                    stop = line[4]
                    strand = line[6]
                    
                    attributes = line[8].split(";")
             
                    for index,inhalt in enumerate(attributes):
                                aufteilung =  str(attributes[index])
                                
                                if aufteilung.startswith("Name="):
                                    if aufteilung.find(",") < 0:
                                        ohnegleich = aufteilung.split("=")[1]
                                        name = ohnegleich
                                    else: 
                                        ohnegleich = aufteilung.split("=")[1]
                                        ohnegleich = "\"" + ohnegleich + "\""
                                        name = ohnegleich
                                if aufteilung.startswith("ID="):
                                     if aufteilung.find(",") < 0:
                                        ohnegleich = aufteilung.split("=")[1]
                                        idname = ohnegleich
                                     else: 
                                        ohnegleich = aufteilung.split("=")[1]
                                        ohnegleich = "\"" + ohnegleich + "\""
                                        idname = ohnegleich
                    
                    zeilenmerkmale = [idname, name, start, stop, strand] 
                    
                    #spciehrt den namen des plasmids und danach alle cds, die darin vorkommen
                    merkmalsammlung.append(zeilenmerkmale)
                
                #wenn plasmid schonmal dran kam, einfach weiter infos in merkmalsammlung speichern bis ein neues kommt
                else:
                    start = line[3]
                    stop = line[4]
                    strand = line[6]
                    
                    attributes = line[8].split(";")
             
                    for index,inhalt in enumerate(attributes):
                                aufteilung =  str(attributes[index])
                                
                                if aufteilung.startswith("Name="):
                                    if aufteilung.find(",") < 0:
                                        ohnegleich = aufteilung.split("=")[1]
                                        name = ohnegleich
                                    else: 
                                        ohnegleich = aufteilung.split("=")[1]
                                        ohnegleich = "\"" + ohnegleich + "\""
                                        name = ohnegleich
                                if aufteilung.startswith("ID="):
                                     if aufteilung.find(",") < 0:
                                        ohnegleich = aufteilung.split("=")[1]
                                        idname = ohnegleich
                                     else: 
                                        ohnegleich = aufteilung.split("=")[1]
                                        ohnegleich = "\"" + ohnegleich + "\""
                                        idname = ohnegleich
                    
                    zeilenmerkmale = [idname, name, start, stop, strand] 
                    merkmalsammlung.append(zeilenmerkmale)
                    
    #letztes plasmid ebenfalls abspeichern, da dies im loop noch nicht ins dict aufgenommen wurde sondern nur in merkmalsammlung            
    dictio[merkmalsammlung[0]] = merkmalsammlung[1:]
    return(dictio)
                    

#COG categorien suchen und der ID zuordnen
def COG_cat(file):
    cog_ID = {}
    cog_counter = {}
    for line in file:
        for element in line:
            #wenn die zeile mit ID beginnt
            if element.startswith('ID'):
                #teile die zeile an ; auf 
                items = element.split(";")
                ID = element.split(";")[0]
                for objekt in items:
                    if objekt.startswith('Dbxref'):
                        COG = objekt.split(",")
                        for cog_part in COG:
                            #guckt ob letzter index auch ein buchstabe ist um Dbxref=COG:C,COG:COG1359 nicht mit einzubeziehen sondern nur: Dbxref=COG:COG0583,COG:K 
                            if cog_part.startswith("COG") and cog_part[-1].isalpha():
                                #key/ value trennung bei :
                                key = cog_part.split(":")[0]
                                value = cog_part.split(":")[1]
                                #wenn key COG ist
                                if key == "COG":
                                    #schreibe in die das neue dict
                                    cog_ID.update({ID:value})
                                # prÃ¼fen ob value nur ein zeichen hat oder mehr
                                if len(value) < 2:
                                    try:
                                        cog_counter[value] += 1    
                                    except KeyError:
                                        cog_counter[value] = 1
                                 # hat value mehr zeichen soll es aufgeteilt und getrennt gezÃ¤hlt werden
                                elif len(value) >= 2:
                                    value2 = list(value)
                                    for buchstabe in value2:
                                        try:
                                            cog_counter[buchstabe] += 1    
                                        except KeyError:
                                            cog_counter[buchstabe] = 1
                                else:
                                    pass
                            # verwertet auch EintrÃ¤ge der Art: Dbxref=COG:C,COG:COG1359     
                            elif cog_part.startswith("Dbxref") and cog_part[-1].isalpha():
                                cog = cog_part.split("=")[1]
                                if cog.startswith("COG"):
                                    key = cog.split(":")[0]
                                    value = cog.split(":")[1]
                                    #wenn key COG ist
                                    if key == "COG":
                                        #schreibe in die das neue dict
                                        cog_ID.update({ID:value})
                                    if len(value) < 2:
                                        try:
                                            cog_counter[value] += 1
                                        except KeyError:
                                            cog_counter[value] = 1
                                          
                                    elif len(value) >= 2:
                                        value2 = list(value)
                                        for buchstabe in value2:
                                            try:
                                                cog_counter[buchstabe] += 1    
                                            except KeyError:
                                                cog_counter[buchstabe] = 1
                                 
                                

    print("\nAnzahl der Gene mit COG-Kategorien:")
    print(cog_counter)
    return cog_ID,cog_counter

##### ACHTUNG #### nachher nochmal anschauen ob ein tab oder zwei, jenachdem ob in Excel oder nicht
# COG Categorien interpretieren und in neue Datei ausgeben; name der neu erstellten datei kann mitgegeben werden
def COG_meaning(dictionary, file_name):
    dateiname = file_name.split(".")[0]
    file = open("COG_ID_meaning_"+dateiname+".csv", 'w')
    file.write("key,"+"COG category," + "category meaning\n"+"\n")
    for key,value in dictionary.items():
        # interpretation der COG categorien fÃ¼r einzelne values
        if len(value) < 2:
            if value == "A":
                file.write(key+","+ value+"," + "RNA processing and modification\n")
            elif value == "B":    
                file.write(key+","+ value+"," + "chromatin structure and dynamics\n")
            elif value == "C":
                file.write(key+","+ value+"," + "energy production and conversation\n")
            elif value == "D":
                file.write(key+","+ value+"," + "cell cycle controll/cell division/chromosome partitioning\n")
            elif value == "E":
                file.write(key+","+ value+"," + "amino acid transport and metabolism\n")
            elif value == "F":
                file.write(key+","+ value+"," + "nucleotid transport and metabolism\n")
            elif value == "G":
                file.write(key+","+ value+"," + "carbohydrate transport and metabolism\n")
            elif value == "H":
                file.write(key+","+ value+"," + "coenzyme transport and metabolism\n")
            elif value == "I":
                file.write(key+","+ value+"," + "lipid transport and metabolism\n")
            elif value == "J":
                file.write(key+","+ value+"," + "translation, ribosomal structure and biogensis\n")
            elif value == "K":
                file.write(key+","+ value+"," + "transkription\n")
            elif value == "L":
                file.write(key+","+ value+"," + "replication/recombination/repair\n")
            elif value == "M":
                file.write(key+","+ value+"," + "cell wall/membrane/envelope biogensis\n")
            elif value == "N":
                file.write(key+","+ value+"," + "cell motility\n")
            elif value == "O":
                file.write(key+","+ value+"," + "posttranslational modification/protein turnover/chaperones\n")
            elif value == "P":
                file.write(key+","+ value+"," + "inorganic ion transport and metabolism\n")
            elif value == "Q":
                file.write(key+","+ value+"," + "secondary metabolites biosynthesis/transport/catabolism\n")
            elif value == "R":
                file.write(key+","+ value+"," + "general function prediction only\n")
            elif value == "S":
                file.write(key+","+ value+"," + "funktion unknown\n")
            elif value == "T":
                file.write(key+","+ value+"," + "signal transduction mechanism\n")
            elif value == "U":
                file.write(key+","+ value+"," + "intercellular trafficking/secretion/vesicular transport\n")
            elif value == "V":
                file.write(key+","+ value+"," + "defense mechanisms\n")
            elif value == "W":
                file.write(key+","+ value+"," + "extracellular structures\n")
            elif value == "X":
                file.write(key+","+ value+"," + "mobilome: prophages/transposons\n")
            elif value == "Y":
                file.write(key+","+ value+"," + "nuclear structure\n")
            elif value == "Z":
                file.write(key+","+ value+"," + "cytoskeleton\n")
                
        # interpretation der COG categorien fÃ¼r doppelte values z.B: OE usw.
        elif len(value) >= 2:
            COG_einzeln =list(value)
            for categorie in COG_einzeln:
                if categorie == "A":
                    file.write(key+","+ value+"," + "RNA processing and modification\n")
                elif categorie == "B":    
                    file.write(key+","+ value+"," + "chromatin structure and dynamics\n")
                elif categorie == "C":
                    file.write(key+","+ value+"," + "energy production and conversation\n")
                elif categorie == "D":
                    file.write(key+","+ value+"," + "cell cycle controll/cell division/chromosome partitioning\n")
                elif categorie == "E":
                    file.write(key+","+ value+"," + "amino acid transport and metabolism\n")
                elif categorie == "F":
                    file.write(key+","+ value+"," + "nucleotid transport and metabolism\n")
                elif categorie == "G":
                    file.write(key+","+ value+"," + "carbohydrate transport and metabolism\n")
                elif categorie == "H":
                    file.write(key+","+ value+"," + "coenzyme transport and metabolism\n")
                elif categorie == "I":
                    file.write(key+","+ value+"," + "lipid transport and metabolism\n")
                elif categorie == "J":
                    file.write(key+","+ value+"," + "translation/ribosomal structure/biogensis\n")
                elif categorie == "K":
                    file.write(key+","+ value+"," + "transkription\n")
                elif categorie == "L":
                    file.write(key+","+ value+"," + "replication/recombination/repair\n")
                elif categorie == "M":
                    file.write(key+","+ value+"," + "cell wall/membrane/envelope biogensis\n")
                elif categorie == "N":
                    file.write(key+","+ value+"," + "cell motility\n")
                elif categorie == "O":
                    file.write(key+","+ value+"," + "posttranslational modification/protein turnover/chaperones\n")
                elif categorie == "P":
                    file.write(key+","+ value+"," + "inorganic ion transport and metabolism\n")
                elif categorie == "Q":
                    file.write(key+","+ value+"," + "secondary metabolites biosynthesis/transport/catabolism\n")
                elif categorie == "R":
                    file.write(key+","+ value+"," + "general function prediction only\n")
                elif categorie == "S":
                    file.write(key+","+ value+"," + "funktion unknown\n")
                elif categorie == "T":
                    file.write(key+","+ value+"," + "signal transduction mechanism\n")
                elif categorie == "U":
                    file.write(key+","+ value+"," + "intercellular trafficking/secretion/vesicular transport\n")
                elif categorie == "V":
                    file.write(key+","+ value+"," + "defense mechanisms\n")
                elif categorie == "W":
                    file.write(key+","+ value+"," + "extracellular structures\n")
                elif categorie == "X":
                    file.write(key+","+ value+"," + "mobilome: prophages/transposons\n")
                elif categorie == "Y":
                    file.write(key+","+ value+"," + "nuclear structure\n")
                elif categorie == "Z":
                    file.write(key+","+ value+"," + "cytoskeleton\n")
            
    file.close()
    #print("\n")
    logging.info("\nCSV mit COG-Kategorien erstellt")

#dna Sequenz in rna umwandeln
def dna_rna(dna_seq):
    valid_dna = "atgcn"
    lower_mRNA ={}
    for key, value in dna_seq.items():
        for nukleotid in value:
            # überprüft ob es sich überhaupt um eine geeignete dna sequenz handelt
            if not nukleotid in valid_dna:
                raise Exception("ERROR: " + str(nukleotid) + " is no valid DNA sequence")
            # wenn es eine geeigente dna sequenz ist wird sie in rna (also ohne t's) umgeschrieben
            else:
                if nukleotid == "t":
                    try:
                        lower_mRNA[key]+= "u"
                    except KeyError:
                        lower_mRNA[key]= "u"
                else:
                    try:
                        lower_mRNA[key]+= nukleotid
                    except KeyError:
                        lower_mRNA[key]=nukleotid
                   
    return lower_mRNA

# Funktion um Gene auf dem - strang umzudrehen und zu komplementieren
def komplementieren(mRNA_seq, start, stop):
    kompl_seq = ""
    # bei einer range von 1:15 werden die basen 1-14 abgebildet; rückwärts genauso - deshalb nicht [15:1:-1] sondern [14:0:-1] damit alle von 14-1 abgebildet werden
    # deshalb von start und von stop jeweils 1 abziehen
    start1 = int(start-1)
    stop1 = int(stop-1)
    for base in mRNA_seq[stop1:start1:-1]:
        # basen komplementär umschreiben und neu in string abspeichern
        if base ==  "t":
            kompl_seq += "a"
        elif base == "a":
            kompl_seq += "t"
        elif base == "g":
            kompl_seq += "c"
        elif base == "c":
            kompl_seq += "g"
    return kompl_seq

# Nukleotidsequenz in AS-Sequenz übersetzen
def translate(mRNA_seq):
    aaseq = ""
    for nukleotid in range(0, len(mRNA_seq), 3):
        # festlegen dass die basen immer in 3er schritten gelesen werden sollen
        triplet = mRNA_seq[nukleotid:nukleotid+3]   
        if (triplet == "att" or triplet == "atc" or triplet == "ata"):
            aaseq += "I"
        elif (triplet == "atg"):
            aaseq += "M"
        elif (triplet in ["act", "acc", "aca", "acg"]):
            aaseq += "T"
        elif (triplet in ["aat", "aac"]):
            aaseq += "N"
        elif (triplet in ["aaa", "aag"]):
            aaseq += "K"
        elif (triplet in ["agt", "agc", "tct", "tcc", "tca", "tcg"]):
            aaseq += "S"   
        elif (triplet in ["agg", "aga", "cgt", "cgc", "cga", "cgg"]):
            aaseq += "R"
        elif (triplet in ["gtt", "gtc", "gta", "gtg"]):   
            aaseq += "V"
        elif (triplet in ["gct", "gcc", "gca", "gcg"]):
            aaseq += "A"
        elif (triplet in ["gat", "gac"]):
            aaseq += "D"
        elif (triplet in ["gaa", "gag"]):
            aaseq += "E"
        elif (triplet in ["ggt", "ggc", "gga", "ggg"]):
            aaseq += "G"
        elif (triplet in ["ttt", "ttc"]):
            aaseq += "F"
        elif (triplet in ["tta", "ttg", "ctt", "ctc", "cta", "ctg"]):
            aaseq += "L"
        elif (triplet in ["tat", "tac"]):
            aaseq += "Y"
        elif (triplet in ["taa", "tag", "tga"]):
            aaseq += "*"
        elif (triplet in ["tgt", "tgc"]):
            aaseq += "C"
        elif (triplet in ["tgg"]):
            aaseq += "W"
        elif (triplet in ["cct", "ccc", "cca", "ccg"]):
            aaseq += "P"
        elif (triplet in ["cat", "cac"]):
            aaseq += "H"
        elif (triplet in ["caa", "cag"]):
            aaseq += "Q"  
    return aaseq


# funktion die alles zusammenführt und auswertet; gibt ein dict aus mit genID als key und AS sequenz als value
def auswerten(gff3, fasta, file_name):
    dateiname = file_name.split(".")[0]
    nukseq = ""
    rna_fasta = fasta
    # Dateien öffnen in die geschrieben werden soll
    filenuk = open("CDS_NukSeq_"+dateiname+".fna", 'w')
    fileaas = open("CDS_AASeq_"+dateiname+".faa", 'w')
    filenuk.write("#Plasmid\t"+"protein-ID\t"+"Nukleotidsequenz\n")
    fileaas.write("#Plasmid\t"+"protein-ID\t"+"Aminosäuresequenz\n")
    # auf keys und values von der in rna umgewandelten fasta zugreifen
    for key_rna_fasta, value_rna_fasta in rna_fasta.items():
        for key_gff3, value_gff3 in gff3.items():
            #werte_gff3[0] = die gen ID; 1 ist start; 2 ist stop 3 ist + oder -
            for werte_gff3 in value_gff3:
                # gen ID abspeichern
                gen_ID = werte_gff3[0]
                start = int(werte_gff3[2])
                stop = int(werte_gff3[3])
                nukseq = ""
                
                for char in range(start,stop):
                    if stop < len(value_rna_fasta):
                        nukseq += value_rna_fasta[char]
                
                parts = [nukseq[i:i+80] for i in range(0, len(nukseq)-1, 80)]
                nukseq = '\n'.join(parts)

                # wenn das gen auf dem "-" strang liegt, muss es komplemetiert werden bevor es translatiert werden kann weil in der fasta ja alles in 5'-3' steht 
                if werte_gff3[4] == "-":
                    #start und stop für die komplementieren funktion abspeichern

                    # sicherstellen dass die sequenz für das gen vom richtigen plasmid stammt
                    if key_gff3 == key_rna_fasta:
                        # sequenzabschnitt komplementieren von start zu stop 
                        komp_seq = komplementieren(value_rna_fasta, start, stop)
                        # sequenzabschnitt translatieren
                        AA_seq = translate(komp_seq)
                        
                        teile = [AA_seq[i:i+80] for i in range(0, len(AA_seq)-1, 80)]
                        AA_seq = '\n'.join(teile)
                        # übersetzte sequenz wieder zusammen mit dem key(nummer des plasmids) und der protein ID abspeichern
                        fileaas.write(">"+key_gff3+"\t" +gen_ID +"\n"+ AA_seq+"\n")
                        filenuk.write(">"+key_gff3+"\t" +gen_ID +"\n"+ nukseq.upper() +"\n")
                # wenn das gen auf dem "+" strang liegt, muss es nicht komplementiert werden sondern kann einfach abgelesen werden
                elif werte_gff3[4] == "+":
                    # wieder start und stop festlegen, aber diesmal für die translations funktion, damit nur der ausschnitt übersetzt wird
                   
                    if key_gff3 == key_rna_fasta:
                        AA_seq_plus = translate(value_rna_fasta[start:stop])
                         
                        teile = [AA_seq_plus[i:i+80] for i in range(0, len(AA_seq_plus)-1, 80)]
                        AA_seq_plus = '\n'.join(teile)
                        # übersetzte sequenz wieder zusammen mit dem key(nummer des plasmids) und der protein ID abspeichern
                        fileaas.write(">"+key_gff3+"\t" +gen_ID +"\n"+ AA_seq_plus+"\n")
                        filenuk.write(">"+key_gff3+"\t" +gen_ID +"\n"+ nukseq.upper() +"\n")
    print("\n")
    logging.info("\nMultiple Fasta mit jeweils Nukleotid- und Aminosäuresequenzen aus der GFF3 erstellt")
    filenuk.close()
    fileaas.close()
    #return file



#Grafik mit einem Dict aus den gewonnenen Daten erstellen
def visualize(daten,dateiname,anhang):
    x = []
    y = []
    xlabel = "Merkmal"
    ylabel = "Anzahl"
    title = anhang
    name = dateiname.split(".")[0] + "_Graph"
    rot=90
    
    for key,value in daten.items():
            
        x.append(key)
        y.append(value)
        
        if len(key) == 3:
            xlabel = "Codon"
        else:
           xlabel = "Merkmal"
           ylabel = "Anzahl" 
           rot=45

    # dataframe = pd.DataFrame(data=daten, index= idx)
    # dataframe.plot(kind = "bar", legend = True)
    plt.figure(figsize=(20,20))
    plt.xticks(rotation=rot)
    plt.xticks(size=15)
    plt.yticks(size=15)
    plt.bar(x,y)
    plt.xlabel(xlabel,fontsize=30)
    plt.ylabel(ylabel, fontsize = 30)
    plt.title(title, fontsize = 30)
    
    plt.savefig(title+"_"+name+'.png')
    #print("\n")
    logging.info("\n"+title+"-Plot erstellt")


def main():
    #PARSER
    parser = argparse.ArgumentParser(description="Programm zur Analyse von gff3 Dateien von Katharina Steeg und Moritz Hobein \n Ausgaben: ")
    parser.add_argument('File', help="Dateiname, der eingelesen werden soll. GFF3-Datei")
    args = parser.parse_args()
    
    
    
    
    
    #Durchführung der Funktionen wie durch argpars eingestellt
    
    Datei = open(args.File,"r")

    Dateiname = os.path.basename(Datei.name)
    meine_Datei,fasta = speichern(Datei)
    merkmale = suchen(meine_Datei)
    
    übersicht(meine_Datei,Dateiname)
    
    COG_ID,anzahl_categories = COG_cat(meine_Datei)
  
    COG_meaning(COG_ID, Dateiname)
    
    gene = genextraktion(meine_Datei)
    
    auswerten(gene, fasta, Dateiname)
    
    codons,startcodons = codoncount(fasta,meine_Datei)
    

    
    visualize(merkmale,Dateiname,"Merkmale")
    visualize(codons,Dateiname,"Codonusage")
    visualize(startcodons,Dateiname,"Startcodons")
    
  
    
    
    



if __name__ == '__main__':
    main()
