# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 20:12:30 2021

@author: moritz hobein & katharina steeg
"""
import os.path
import matplotlib.pyplot as plt
import argparse
import logging

# Ursprungsdatei für die Weiterverarbeitung vorbereiten
def speichern(file, fastastate):
    
    #Überprüfen, ob die Datei sich selbst als gff3-File bezeichnet
    if not file.readline().startswith("##gff-version 3"):
        logging.error("Es wurde scheinbar kein gff3-File eingelesen")
        
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
    
    #für fasta-Aufbereitung; erfolgt nur hier, wenn keine zweite fasta übergeben wurde. Diese Information wird über fastastate vermittelt
    if not fastastate:
        file.seek(0)
        plasmide = separateFasta(file)
        return liste,plasmide
    
    else:
        return liste

#fasta-Aufbereitung, wenn eine separate fasta-Datei eingelesen wurde. Sollte dennoch eine an die Sequenz an die gff3-Datei angehängt sein, wird dennoch diese hier genommen    
def separateFasta(fasta):
    plasmide = {}
    anzahl = []
    zähler = 0
    seq = ""
    ohnev = ""

    for line in fasta:
        
        #guckt, ob eine neue separate Sequenz anfängt
        if line.endswith("\n") and line.startswith(">"):
            zeile = line.strip()
            
            #String mit Sequenz wird zusammen mit Sequenznamen ins dict eingetragen
            plasmide[ohnev] = seq.lower()
            
            #Stringreset für nächte Sequenz
            seq = ""
            
            
            ohnev = zeile.split(">")[1]
            zähler += 1 
            anzahl.append(line)
        
        #speichert alle Zeilen, die nicht der Sequenzname ist, ab
        if line.endswith("\n") and zähler > 0:
            if line.startswith(">"):
                pass
            else:
                zeile  = line.strip()
                
                # Zähler zählt mit, wie viele separate Sequenzen da sind. Wenn eine neue beginnt, wird der string seq resettet, um die nächste Sequenz abzuspeichern
                # seq wird nur erweitert, bis eine neue seq beginnt (der Zähler ist dann nicht mehr gleich Anzahl, weil dieses um einen Eintrag erweitert wurde)
                if zähler == len(anzahl):
                     
                    #speichert jede Zeile aneinandergehängt in einen String
                    seq += zeile
    
    #weil am Ende des files zwar ohnev und seq mit den Daten der letzten Sequenz gefüllt sind, aber der loop vor dem letzten Eintrag ins dict aufhört, kommt der Eintrag hier
    plasmide[ohnev] = seq.lower()
    
    #entfernen des leeren Eintrags am Anfang, der durch die Funktionsstruktur bedingt war
    plasmide.pop("")  
    
    
    return plasmide


#Suche nach Merkmalen und ihrer Anzahl
def suchen(file, cmd):
    zähler = {}
    for line in file:
        # alle Zeilen überspringen, die nicht im Tabellenformat der gff3-Datei sind mit 9 verschiedenen Fields
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
            
            #Zählen der Merkmale mit COG category
            if "Dbxref" in aufteilung:
                    try:
                        aufteilung.find("COG:")
                        zähler["Merkmale mit COG category"] +=1
                    except KeyError:
                        zähler["Merkmale mit COG category"] = 1
                    except :
                        pass
    if cmd:           
        print("\nannotierte Sequenzabschnitte samt Anzahl:")          
        print(sorted(zähler.items(), key = lambda kv: kv[1], reverse=True),'\n')              
    
    return zähler


# def ListToString(liste):
#     ausgabe = ""
#     for element,index in enumerate(liste):
#         ausgabe += liste[element]+"\t"
#     print(ausgabe)
#     return ausgabe


#Suchen nach der Häufigkeit des Vorkommens eigener Terme, die über argparse mitgegeben wurden 
def manuellSuchen(file, gesuchtes, ignorieren="hierkönnteihrewerbungstehen!!11!1!"):
    zähler = {}
    
    #zum Ignorieren von Zeilen beim Zählen, wenn der in ignorieren abgespeicherte Term auch in der Zeile vorkommt, in der gesuchtes gefunden wurde
    exception = False
    
    #für die Erzeugung eines verkürzten Teils der Annotationszeilen, falls diese exportiert werden sollen
    zeilen = []
    
    #alle Spalten bis auf "Attribute", da hier die Ausgabe sehr unübersichtlich und ich denke auch nicht zielführend ist. Kann man natürlich auf Wunsch auch ändern. 4 und 5 sind auch nicht sonderlich hilfreich.
    fields = "123456789"

    #sollte man ein - oder ein + suchen, wird davon ausgegangen, dass es sich um eine Suche nach der Anzahl von + und - Strängen handelt (- als Bindestrich kommt sicherlich häufig vor, ist aber denke ich eine recht irrelevante Information und wird deswegen hiermit ausgesiebt)
    if gesuchtes == "+" or gesuchtes == "-":
        strang = True
    else:
        strang = False
        
    #gibt man eine der Zahlen aus fields ein, werden die Inhalte der betreffenden Spalte gezählt
    if str(gesuchtes) in fields: 
        field = True
        
        #Name der Spalte wird hier für die Ausgabe gespeichert
        if gesuchtes == "1":
            fieldname = "SeqID"
            
        elif gesuchtes == "2":
            fieldname = "Source"
        
        elif gesuchtes == "3":
            fieldname = "Type"
        
        elif gesuchtes == "4":
            fieldname = "Start"
            
        elif gesuchtes == "5":
            fieldname = "End"
        
        elif gesuchtes == "6":
            fieldname = "Score"
        
        elif gesuchtes == "7":
            fieldname = "Strand"
        
        elif gesuchtes == "8":
            fieldname = "Phase"
        
        elif ignorieren == "9":
            fieldname = "Attributes"
            
    else:
        field = False
    
    #gleiches wie vorher nur für die Terme, die ausgeschlossen werden sollen
    if ignorieren == "+" or ignorieren == "-":
        exceptstrang = True
    else:
        exceptstrang = False
        
    if ignorieren in fields:
        exceptfield = True
        #Name der Spalte wird hier für die Ausgabe gespeichert
        if ignorieren == "1":
            exfieldname = "SeqID"
            
        elif ignorieren == "2":
            exfieldname = "Source"
        
        elif ignorieren == "3":
            exfieldname = "Type"
        
        elif ignorieren == "4":
            exfieldname = "Start"
            
        elif ignorieren == "5":
            exfieldname = "End"
        
        elif ignorieren == "6":
            exfieldname = "Score"
        
        elif ignorieren == "7":
            exfieldname = "Strand"
        
        elif ignorieren == "8":
            exfieldname = "Phase"
            
        elif ignorieren == "9":
            exfieldname = "Attributes"
    else:
        exceptfield = False
            
    #Durchsuchen der einzelnen Zeilen
    for line in file:
        # alle Zeilen überspringe, die nicht im Tabellenformat der gff3-Datei sind mit 9 verschiedenen Fields
        if len(line) == 9:
            zeilenzähler = 0
            
            #wenn + oder - eingegeben wurde, nur im für den Strang betreffenden Field nachschauen
            if strang:
                
                if exceptstrang:
                    logging.error("\nSuche ergibt keinen Sinn, versuche etwas anderes (Strang-ohne-Strang-Suche)\n")
                    break
                
                #Warunungen, wenn sich die Eingaben gegenseitig ausschließen
                elif exceptfield:
                    if ignorieren == "6":
                       logging.error("\nSuche ergibt keinen Sinn, versuche etwas anderes (Strang-ohne-Strang-Suche)\n")
                       break
                    else:
                        logging.error("\n Suche ergibt keinen Sinn, versuche etwas anderes (Strang kann nicht ohne eine andere Spalte existieren)")
                        break
                
                #Suchen der Anzahl und Zeilen des gewünschten Stranges, der ein Merkmal nicht vorweist
                else:
                    if line[6] == gesuchtes:
                        for element in line:
                        
                            if element == ignorieren:
                                exception = True
                            #letzte Spalte mit den Attributen wird gesplittet, sodass die Attribute einzeln in einer Liste stehen
                        attributes = line[8].split(";")
                        
                        #dann nach dem Term im letzten Field suchen
                        for element in attributes:
                            
                        
                            if element.find(ignorieren) > 0:
                                exception = True
                        
                        if not exception:
                            #wenn etwas gefunden wird, wird die Zeile als ganzes abgespeichert für eine eventuelle spätere Ausgabe 
                            zeilen.append(line)
                            try:
                               zähler[gesuchtes] += 1
                            except KeyError:
                               zähler[gesuchtes] = 1
                exception = False
                
            #eine einzelne Spalte auszählen
            elif field:
                if exceptfield:
                    
                   #Warnung ausgeben, wenn sich die Eingaben gegenseitig ausschließen
                   if gesuchtes == ignorieren:
                       logging.error("\nSuche ergibt keinen Sinn, versuche etwas anderes (Suche und Ausschluss im gleichen Field)\n")
                       break
                   else:
                       logging.error("\nSuche ergibt keinen Sinn, versuche etwas anderes (Spalte kann nicht ohne eine andere Spalte existieren)")
                       break
                
                #Die Inhalte einer Spalte ohne eine bestimmte Strangart wird mit Häufigkeiten des Vorkommens und der Zeile abgespeichert
                elif exceptstrang:
                    if line[6] == ignorieren:
                        exception = True
                    else:
                       #alle Spalten außer "Attribute"
                        if gesuchtes != "9":
                        
                            if not exception:
                                    #wenn etwas gefunden wird, wird die Zeile als ganzes abgespeichert für eine eventuelle spätere Ausgabe 
                                    zeilen.append(line)
                                    try:
                                       zähler[line[int(gesuchtes)-1]] += 1
                                    except KeyError:
                                       zähler[line[int(gesuchtes)-1]] = 1 
                
                #Die Inhalte iner Spalte ohne einen bestimmten Term werden gespeichert
                else:
                    #alle Spalten außer "Attribute"
                    if gesuchtes != "9":
                        for element in line:
                        
                            if element == ignorieren:
                                exception = True
                            
                            #letzte Spalte mit den Attributen wird gesplittet, sodass die Attribute einzeln in einer Liste stehen
                        attributes = line[8].split(";")
                        
                        #dann nach dem Term im letzten Field suchen
                        for element in attributes:

                            if element.find(ignorieren) > 0:
                                exception = True
                            
                        if not exception:
                                #wenn etwas gefunden wird, wird die Zeile als ganzes abgespeichert für eine eventuelle spätere Ausgabe 
                                zeilen.append(line)
                                try:
                                   zähler[line[int(gesuchtes)-1]] += 1
                                except KeyError:
                                   zähler[line[int(gesuchtes)-1]] = 1
                exception = False
            
            #alle anderen Suchanfragen
            else:
                
                #Zeilen, in denen der Term ohne die gewünschte Strangart vorkommt, werden gespeichert
                if exceptstrang:
                  if line[6] == ignorieren:
                      exception == True
                 
                  else:
                      for element in line:
                    
                        #erst Nachschauen, ob eines der ersten 8 Fields genau wie der Term heißt
                        if element == gesuchtes:
                             zeilenzähler +=1
                        
                        #letzte Spalte mit den Attributen wird gesplittet, sodass die Attribute einzeln in einer Liste stehen
                        attributes = line[8].split(";")
                        
                        #dann nach dem Term im letzten Field suchen
                        for element in attributes:
                            
                        
                            if element.find(gesuchtes) > 0:
            
                                    zeilenzähler += 1
                                    
                #Zeilen, in denen der gesuchte Term außerhalb der gewünschten Spalte vorkommt
                elif exceptfield:
                    
                    #für Spalten 1-8
                    if line[int(ignorieren)-1] == gesuchtes:
                        exception = True
                    
                    #für Spalte 9
                    elif ignorieren == "9" and line[int(ignorieren)-1].find(gesuchtes) > 0:
                        exception = True
                   
                    else:
                        for element in line:
                    
                            #erst Nachschauen, ob eines der ersten 8 Fields genau wie der Term heißt
                            if element == gesuchtes:
                                 zeilenzähler +=1
                            
                            #letzte Spalte mit den Attributen wird gesplittet, sodass die Attribute einzeln in einer Liste stehen
                            attributes = line[8].split(";")
                            
                            #dann nach dem Term im letzten Field suchen
                            for element in attributes:
                                
                            
                                if element.find(gesuchtes) > 0:
                
                                        zeilenzähler += 1
                                        
                              
                else:
                    
                    for element in line:
                        
                        #erst Nachschauen, ob eines der ersten 8 Fields genau wie der Term heißt
                        if element == gesuchtes:
                             zeilenzähler +=1
                        
                        if element == ignorieren:
                            exception = True
                        #letzte Spalte mit den Attributen wird gesplittet, sodass die Attribute einzeln in einer Liste stehen
                        attributes = line[8].split(";")
                        
                        #dann nach dem Term im letzten Field suchen
                        for element in attributes:
                            
                        
                            if element.find(gesuchtes) > 0:
            
                                    zeilenzähler += 1
                                    
                            try:
                                element.index(ignorieren)
                                exception = True
                            except:
                                pass

                #Doppelzählungen innerhalb einer Zeile werden nur als eine einzige Sichtung gewertet (kommt ja häufig vor, dass z.B. name und product gleich sind, aber das ist uns ja wahrscheinlich egal)
                if zeilenzähler > 0 and not exception:
                    zeilen.append(line)
                        #Anlegen eines dict mit allen Merkmalen und deren Häufigkeit in der Datei
                    try:
                        zähler[gesuchtes] += 1
                    except KeyError:
                        zähler[gesuchtes] = 1
                exception = False
    
    #Printausgabe wird angepasst, je nach dem ob ein Strang oder ein anderer Term gesucht oder ignoriert wurde
    
    if ignorieren != "hierkönnteihrewerbungstehen!!11!1!":
        if strang:
                print("\nAnzahl der Merkmale auf dem "+gesuchtes+ " Strang ohne "+ignorieren+":")
        elif field:
            if exceptfield:
                print("\nInhalte samt Anzahl in Spalte "+gesuchtes+" ("+fieldname+") "+"ohne Spalte "+ignorieren+" ("+exfieldname+"):")
            elif exceptstrang:
                print("\nInhalte samt Anzahl in Spalte "+gesuchtes+" ("+fieldname+") "+"ohne "+ignorieren+" Strang:")
            else:
                print("\nInhalte samt Anzahl in Spalte "+gesuchtes+" ("+fieldname+") "+"ohne "+ignorieren+":")
        else:
            if exceptfield:
               print("\nAnzahl Sequenzabschnitte, in deren Annotation der gesuchte Term "+gesuchtes+" außerhalb von Spalte "+ignorieren+" ("+exfieldname+") vorkommt:")
            elif exceptstrang:
                print("\nAnzahl Sequenzabschnitte, in deren Annotation der gesuchte Term "+gesuchtes+" außerhalb von "+ignorieren+" Strängen vorkommt:")
            else:
                print("\nAnzahl Sequenzabschnitte, in deren Annotation der gesuchte Term "+gesuchtes+" ohne "+ignorieren+" vorkommt:")
    else:
        if strang:
            print("\nAnzahl der Merkmale auf dem "+gesuchtes+ " Strang:")
        elif field:
            print("\nInhalte samt Anzahl in Spalte "+gesuchtes+" ("+fieldname+"):")
        else:
            print("\nAnzahl Sequenzabschnitte, in deren Annotation der gesuchte Term "+gesuchtes+" vorkommt:")          
    
    print(sorted(zähler.items(), key = lambda kv: kv[1], reverse=True)) 
    print("\n")             
    
    return zähler,zeilen

#Erstellen einer Übersichts-csv mit den wesentlichen Informationen der Merkmale
def übersicht(file, filename):
    
    #Vorbereitung für den Namen der neu erstellten Datei
    ohneendung = filename.split(".")[0]
    outfile = open("Übersicht_"+ ohneendung + ".csv","w")
    
    #Kopfzeile der Datei
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
           
          #Erzeugen der Zeile mit den ganzen Basisinformationen eines Merkmals
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
    logging.info("\nÜbersicht als CSV erstellt\n")

#für die manuelle Suche mit dem Ignorieren von Ergebnissen mit bestimmten Inhalten
def suchenOhne(liste, file):
    zähler,zeilen = manuellSuchen(file, liste[0], liste[1])
   
    return zähler,zeilen 


#Zählt die Anzahl vorkommender Codons in den CDS            
def codoncount(seq, gff3, cmd):
    
    anzahl = {}
    anteil = {}
    insgesamt = 0
    startinsg = 0
    startcodons = {}
    startanteil = {}
    
    for key,value in seq.items():
        
       for element in gff3:
           
           #Speichern der Grenzparameter der CDS; wird nur ausgewertet auf dem richtigen Plasmid durch Abgleichen der Keys
           if len(element) == 9 and key == element[0] and element[2] == "CDS": 
               start = int(element[3])
               stop = int(element[4])
               
               #für den Fall, dass die Grenzparameter der CDS die Sequenzlänge überschreiten, wird das vermerkt und die CDS übersprungen
               if stop > len(value):
                   
                   #Beschaffen der ID der zu langen CDS
                   attributes = element[8].split(";")
                        
                   for index,inhalt in enumerate(attributes):
                       aufteilung =  str(attributes[index])
                
                       if aufteilung.startswith("ID="):
                          ohnegleich = aufteilung.split("=")[1]
                          name = ohnegleich
                          #Informieren des Users über die Überlänge der CDS
                   logging.warning("\nannotierte CDS mit ID "+ name +" überschreitet Sequenzlänge und wird deswegen nicht ausgewertet\n")
               else:
                   
                   #Suchen Codons auf dem Sequenzabschnitt, der durch Start und Stop der CDS angegeben wurde
                   if element[6] == "+":
                           
                       #Schleife, die ein einzelnes Codon pro Durchlauf extrahiert
                       for char in range(start-1,stop-3,3):
                           
                               if char == start-1:
                                   
                                   #Zählen des ersten Codons der CDS als Startcodon
                                   starttriplet = value[char] + value[char+1] + value[char+2]
                                   startcodon = starttriplet.replace("t","u")
                                   try:
                                       startcodons[startcodon] += 1
                                   except KeyError:
                                       startcodons[startcodon] = 1
                            
                               #Zusammensetzen des Codons
                               triplet = value[char] + value[char+1] + value[char+2]
                               
                               #Umwandlung von DNA in RNA
                               codon = triplet.replace("t","u")
                               
                               #Hochzählen des ausgelesenen Codons
                               try:
                                   anzahl[codon] += 1
                               except KeyError:
                                   anzahl[codon] = 1
                   
                   #für CDS auf dem Komplementärstrang, das gleiche noch einmal
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
    
    #Zum Ausrechnen der Prozentualen Anteile                
    for key, value in anzahl.items():
        insgesamt += value
        
    for key, value in startcodons.items():
        startinsg += value
    
    for key, value in anzahl.items():
        anteil[key] = round(value/insgesamt*100,2) 
        
    for key,value in startcodons.items():
        startanteil[key] = round(value/startinsg*100,2)
    
    if cmd:
        print("Codonusage:")
        print(sorted(anzahl.items(), key = lambda kv: kv[1], reverse=True))
        print("\nProzentualer Anteil der genutzten Codons:")
        print(sorted(anteil.items(), key = lambda kv: kv[1], reverse=True))
        print("\nAnzahl der genutzten Startcodons")
        print(sorted(startcodons.items(), key = lambda kv: kv[1], reverse=True))
        print("\nAnteil der genutzten Startcodons von der Gesamtheit:") 
        print(sorted(startanteil.items(), key = lambda kv: kv[1], reverse=True))
    
    return anzahl,startcodons   

#Sucht die Informationen zu Weiterverarbeitung der einzelnen Gene heraus und speichere sie in einem Dictionary
def genextraktion(file):
    
    plasmide = []
    merkmalsammlung = []
    dictio = {}
    
    
    for line in file:
        
        #wenn das Format erkannt wird, in dem CDS annotiert wurden
        if len(line) == 9:
            if line[2] == "CDS":
                
                
                #wenn ein neues Plasmid/eine neue Sequenz gelesen wird
                if line[0] not in plasmide:
                    
                    #letztes Plasmid im dict abspeichern mit Informationen der einzelnen CDS aus merkmalsammlung
                    try:
                        dictio[merkmalsammlung[0]] = merkmalsammlung[1:]
                        merkmalsammlung.clear()
                    except:
                        pass
                    
                    
                    #neues Plasmid merken, anfangen mit dem Abspeichern der Zeileninformationen für erste CDS
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
                    
                    #speichert den Namen des Plasmids und danach alle CDS, die darin vorkommen
                    merkmalsammlung.append(zeilenmerkmale)
                
                #wenn das Plasmid schon einmal dran kam, einfach weiter Infos in Merkmalsammlung speichern bis ein neues kommt
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
                    
    #letztes Plasmid ebenfalls abspeichern, da dies im loop noch nicht ins dict aufgenommen wurde sondern nur in merkmalsammlung            
    
    try:
        dictio[merkmalsammlung[0]] = merkmalsammlung[1:]
    except:
        logging.error('Es konnten keine Gene extrahiert werden, da scheinbar keine Annotationen vorhanden sind.')
    return(dictio)
                    

#COG-Kategorien suchen und der ID zuordnen
def COG_cat(file,cmd):
    cog_ID = {}
    cog_counter = {}
    for line in file:
        for element in line:
            #wenn die Zeile mit ID beginnt
            if element.startswith('ID'):
                #teile die Zeile an ; auf 
                items = element.split(";")
                ID = element.split(";")[0]
                for objekt in items:
                    # wenn die objekt in items mit "Dbxref" beginnt
                    if objekt.startswith('Dbxref'):
                        #teile das objekt an "," auf
                        COG = objekt.split(",")
                        for cog_part in COG:
                            #guckt ob letzter Index auch ein Buchstabe ist um z.B. Dbxref=COG:C,COG:COG1359 nicht mit einzubeziehen sondern nur: Dbxref=COG:COG0583,COG:K 
                            if cog_part.startswith("COG") and cog_part[-1].isalpha():
                                #key/ value Trennung bei :
                                key = cog_part.split(":")[0]
                                value = cog_part.split(":")[1]
                                #wenn key COG ist
                                if key == "COG":
                                    #schreibe ID und value in das dict
                                    cog_ID.update({ID:value})
                                # prüfen, ob die COG categorie nur ein Zeichen hat oder mehr
                                if len(value) < 2:
                                    try:
                                        cog_counter[value] += 1    
                                    except KeyError:
                                        cog_counter[value] = 1
                                 # hat die COG categorie mehr Zeichen soll es aufgeteilt und getrennt gezählt werden
                                elif len(value) >= 2:
                                    value2 = list(value)
                                    for buchstabe in value2:
                                        try:
                                            cog_counter[buchstabe] += 1    
                                        except KeyError:
                                            cog_counter[buchstabe] = 1
                                else:
                                    pass
                            # verwertet auch Einträge der Art: Dbxref=COG:C,COG:COG1359 die oben rausgefiltert wurden indem es überprüft ob das letzte Zeichen ein Buchstabe ist     
                            elif cog_part.startswith("Dbxref") and cog_part[-1].isalpha():
                                cog = cog_part.split("=")[1]
                                if cog.startswith("COG"):
                                    key = cog.split(":")[0]
                                    value = cog.split(":")[1]
                                    #wenn key COG ist
                                    if key == "COG":
                                        #schreibe ID und value in das dict
                                        cog_ID.update({ID:value})
                                    if len(value) < 2:
                                        try:
                                            cog_counter[value] += 1
                                        except KeyError:
                                            cog_counter[value] = 1
                                    # wieder Aufteilung bei COG categorien mit zwei Buchstaben oder mehr      
                                    elif len(value) >= 2:
                                        value2 = list(value)
                                        for buchstabe in value2:
                                            try:
                                                cog_counter[buchstabe] += 1    
                                            except KeyError:
                                                cog_counter[buchstabe] = 1
                                 
                                
    if cmd:                            
        print("\nAnzahl der Gene mit COG-Kategorien:")
        print(sorted(cog_counter.items(), key = lambda kv: kv[1], reverse=True),"\n")
        
    return cog_ID,cog_counter


# COG-Kategorien interpretieren und in neue Datei ausgeben; Name der neu erstellten Datei kann mitgegeben werden
def COG_meaning(dictionary, file_name):
    dateiname = file_name.split(".")[0]
    file = open("COG_ID_meaning_"+dateiname+".csv", 'w')
    file.write("key,"+"COG category," + "category meaning\n"+"\n")
    for key,value in dictionary.items():
        # Interpretation der COG-Kategorien fuer einzelne values
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
                
        # Interpretation der COG-Kategorien für doppelte values z.B: OE usw.
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
    logging.info("\nCSV mit COG-Kategorien erstellt\n")


# Funktion um Gene auf dem - Strang umzudrehen und zu komplementieren
def komplementieren(mRNA_seq, start, stop):
    kompl_seq = ""
    # bei einer range von 1:15 werden die Basen 1-14 abgebildet; rückwärts genauso - deshalb nicht [15:1:-1] sondern [14:0:-1] damit alle von 14-1 abgebildet werden
    # deshalb von start und von stop jeweils 1 abziehen
    start1 = int(start-1)
    stop1 = int(stop-1)
    for base in mRNA_seq[stop1:start1:-1]:
        # Basen revers komplementär umschreiben und neu in string abspeichern
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


# Funktion, die alles zusammenführt und auswertet; gibt eine fasta aus mit Protein-ID; Proteinname; Plasmid und Nukleotid- bzw. Aminosäuresequenz
def auswerten(gff3, fasta, file_name):
    dateiname = file_name.split(".")[0]
    nukseq = ""
    rna_fasta = fasta
    #öffnen der finalen Textdateien
    filenuk = open("CDS_NukSeq_"+dateiname+".fna", 'w')
    fileaas = open("CDS_AASeq_"+dateiname+".faa", 'w')
    filenuk.write("#Protein-ID\tProteinname\tPlasmid\tNukleotidsequenz\n")
    fileaas.write("#Protein-ID\t\tProteinname\tPlasmid\tAminosäuresequenz\n")
    # auf keys und values aus dem genextraktions dictionary zugreifen
    for key_rna_fasta, value_rna_fasta in rna_fasta.items():
        for key_gff3, value_gff3 in gff3.items():
            #werte_gff3[0] = die gen ID; 1 ist start; 2 ist stop 3 ist + oder -
            for werte_gff3 in value_gff3:
                # gen ID abspeichern
                gen_ID = werte_gff3[0]
                start = int(werte_gff3[2])-1
                stop = int(werte_gff3[3])
                nukseq = ""
                
                for char in range(start,stop):
                    if stop < len(value_rna_fasta):
                        nukseq += value_rna_fasta[char]
                
                parts = [nukseq[i:i+80] for i in range(0, len(nukseq)-1, 80)]
                nukseq = '\n'.join(parts)

                # wenn das Gen auf dem "-" Strang liegt, muss es komplemetiert werden bevor es translatiert werden kann weil in der fasta ja alles in 5'-3' steht 
                if werte_gff3[4] == "-":
                    # sicherstellen dass die Sequenz für das Gen vom richtigen Plasmid stammt
                    if key_gff3 == key_rna_fasta:
                        # Sequenzabschnitt komplementieren von Start zu Stop 
                        komp_seq = komplementieren(value_rna_fasta, start, stop)
                        # Sequenzabschnitt translatieren
                        AA_seq = translate(komp_seq)
                        teile = [AA_seq[i:i+80] for i in range(0, len(AA_seq)-1, 80)]
                        AA_seq = '\n'.join(teile)
                        # Übersetzte Sequenz wieder zusammen mit dem key(nummer des plasmids) und der Protein ID abspeichern
                        fileaas.write(">" +gen_ID +"\t"+werte_gff3[1]+"\t"+key_gff3+"\n"+ AA_seq+"\n")
                        filenuk.write(">"+gen_ID +"\t"+werte_gff3[1]+"\t"+key_gff3+"\n"+ nukseq.upper() +"\n")
                # wenn das Gen auf dem "+" Strang liegt, muss es nicht komplementiert werden sondern kann einfach abgelesen werden
                elif werte_gff3[4] == "+":
                    # wieder Start und Stop festlegen, aber diesmal für die translations Funktion, damit nur der Ausschnitt übersetzt wird
                    if key_gff3 == key_rna_fasta:
                        AA_seq_plus = translate(value_rna_fasta[start:stop])
                        teile = [AA_seq_plus[i:i+80] for i in range(0, len(AA_seq_plus)-1, 80)]
                        AA_seq_plus = '\n'.join(teile)
                        # übersetzte Sequenz wieder zusammen mit dem key(nummer des plasmids) und der protein ID abspeichern
                        fileaas.write(">" +gen_ID+"\t"+werte_gff3[1]+"\t"+key_gff3+ "\n"+ AA_seq_plus+"\n")
                        filenuk.write(">"+gen_ID +"\t"+werte_gff3[1]+"\t"+key_gff3+"\n"+ nukseq.upper() +"\n")
    
    logging.info("\nMultiple Fasta mit jeweils Nukleotid- und Aminosäuresequenzen aus der GFF3 erstellt\n")
    filenuk.close()
    fileaas.close()
    #return file

#class zum Erstellen von Graphen
class Graph:
    
    def __init__(self,daten,dateiname,title, rot, figlength, figwidth, titlesize, labelsize, ticksize, textsize):
        self.x = []
        self.y = []
        self.daten = daten
        self. title = title
        self.dateiname = dateiname
        self.name = self.dateiname.split(".")[0] + "_Graph"
        self.ylabel = 'Anzahl'
        self.rot = rot
        self.figlength = figlength
        self.figwidth = figwidth
        self.titlesize = titlesize
        self.textsize = textsize
        self.labelsize = labelsize
        self.ticksize = ticksize
        
        #mitgegebenes Dictionary wird für jede Achse in eine Liste gespeichert.
        for key,value in self.daten.items():
            
            self.x.append(key)
            self.y.append(value)
    
    #bar chart für die Anzahl der Codons und Merkmale
    def bar(self):
        
        if self.title == "Codons":
            xlabel = 'Codons'
        else:
            xlabel = 'Merkmale'
        
        plt.figure(figsize=(self.figlength, self.figwidth))
        plt.xticks(rotation=self.rot)
        plt.xticks(size=self.ticksize)
        plt.yticks(size=self.ticksize)
        plt.bar(self.x, self.y)
        plt.xlabel(xlabel, fontsize=self.labelsize)
        plt.ylabel(self.ylabel, fontsize = self.labelsize)
        plt.title(self.title, fontsize = self.titlesize)
        plt.savefig(self.title+"_"+self.name+'.png')
        logging.info("\n"+self.title+"-Plot erstellt\n")
    
    #pie chart, bietet sich an für die Anteile der Startcodons o.ä.
    def pie(self):
        
        plt.figure(figsize=(self.figlength, self.figwidth))
        plt.pie(self.y ,labels=self.x, autopct='%1.1f%%', textprops={'fontsize': self.textsize})
        plt.title(self.title, fontsize = self.titlesize)
        plt.savefig(self.title+"_"+self.name+'.png')
        logging.info("\n"+self.title+"-Plot erstellt\n")


def main():
    Ausgabeoptionen = ["graph", "cmd", "fasta", "csv"]
    run = True
    ausgabecmd = False
    separatF = False
    
    #PARSER
    parser = argparse.ArgumentParser(description="Programm zur Analyse von gff3-Dateien von Katharina Steeg und Moritz Hobein.\nAusgaben: Anzahl annotierter Merkmale, Codon- und Startcodonusage, Merkmale mit COG-Kategorien sowie CSVs dieser und einer heruntergebrochenen Übersicht der Merkmale, Grafiken.")
    parser.add_argument('File', help="Dateiname, der eingelesen werden soll. gff33-Datei")
    parser.add_argument('-f', '--fasta', dest='sepFasta', default={}, help='Wenn in der gff3-Datei keine FASTA-Sequenz angefügt ist, muss diese separat eingelesen werden.')
    parser.add_argument('-a', '--ausgabe', dest='Ausgabe', default=Ausgabeoptionen, nargs='*', help ="Welche Informationen sollen aus der Datei extrahiert werden? Folgende Optionen bestehen:\n"+str(Ausgabeoptionen))
    parser.add_argument('-m', '--manuell', dest='Manuell', default=[], nargs='*', help="Manuelles Suchen in der Datei nach gewünschten Merkmalen, Anzahl wird ausgegeben\nBesonderheiten: +/- sucht nach Anzahl der jeweiligen Stränge; 1-8 sucht nacht den Inhalten des jeweiligen Fields." )
    parser.add_argument('-l', '--logging', dest='Log', default=20, help='Logging Level setzen. Default = 20. Mögliche Level sind x*10 mit 0 <= x <= 5. Höhere Level zeigen nur die jeweils wichtigeren Informationen.')
    parser.add_argument('-o', '--ohne', dest='Ohne', default=[], nargs=2, help='Gibt Anzahl der Zeilen aus, in denen das erste Merkmale vorkommet, das zweite jedoch gleichzeitig nicht. Es gelten die gleichen besonderheiten wie beim manuellen Suchen')
    parser.add_argument('-e', '--export', dest ='Export', default=False, action='store_true', help='Exportiert die Annotationszeilen, in denen manuell gesuchte Merkmale vorkommen, als CSV.')
    args = parser.parse_args()
    
    
    #Überprüfen, ob gesetzte Optionen Sinn ergeben. Wenn nicht, wird das Programm nicht ausgeführt
    for elemente in range(0,len(args.Ausgabe)):
        
        if args.Ausgabe[elemente] in Ausgabeoptionen:
            pass
        else: 
            run = False
            logging.error("Ungültige Ausgabeoption gesetzt, bitte korrigieren. Eine Liste gültiger Optionen befindet sich auf der Hilfeseite")
    
    logging.basicConfig(level=int(args.Log))
    
    if run:
        #Durchführung der Funktionen wie durch argpars eingestellt
        
        Datei = open(args.File,"r")
    
        Dateiname = os.path.basename(Datei.name)
        
        #Überprüfen, ob eine separate Fasta mitgegeben wurde
        if len(args.sepFasta) == 0:
            meine_Datei,fasta = speichern(Datei,separatF)
        
        else:
            separatF = True
            FastaFile = open(args.File,"r")
            meine_Datei = speichern(Datei,separatF)
            fasta = separateFasta(FastaFile)
        
        
        gene = genextraktion(meine_Datei)

        #hier werden die einzelnen Auswertungen nur dann erzeugt, wenn das wie über argparse mitgeteilt erwünscht ist (keine Angabe -> alles ausgeben)
        if "cmd" in args.Ausgabe:
            ausgabecmd = True
    
        merkmale = suchen(meine_Datei,ausgabecmd)    
        codons,startcodons = codoncount(fasta,meine_Datei,ausgabecmd)
        COG_ID,anzahl_categories = COG_cat(meine_Datei,ausgabecmd)
        
        
        if not args.Export:
            for element in args.Manuell:
                manuellSuchen(meine_Datei, element)
        
        if "csv" in args.Ausgabe:
            
            übersicht(meine_Datei,Dateiname)
            COG_meaning(COG_ID, Dateiname)
            
            if args.Export:
                
                if args.Manuell != []:
                    for element in args.Manuell:
                        Anzahl,kurzeAnnotation = manuellSuchen(meine_Datei, element)
                        exportname = element+"_"+Dateiname
                        übersicht(kurzeAnnotation,exportname)
                        
                if args.Ohne != []:
                    exportname = args.Ohne[0]+"_ohne_"+args.Ohne[1]+"_"+Dateiname
                    zahlen, annotationszeilen = suchenOhne(args.Ohne, meine_Datei)
                    übersicht(annotationszeilen, exportname)
        
        if "fasta" in args.Ausgabe:
            
            auswerten(gene, fasta, Dateiname)
            
        if "graph" in args.Ausgabe:
            
            GMerkmale = Graph(merkmale,Dateiname,"Merkmale",50,40,40,60,40,20,30)
            GCodons = Graph(codons,Dateiname,"Codonusage",90,40,40,60,40,20,30)
            GStartcodons = Graph(startcodons,Dateiname,"Startcodons",25,40,40,60,40,20,30)
             
            GMerkmale.bar()
            GCodons.bar()
            GStartcodons.pie()

  
if __name__ == '__main__':
    main()