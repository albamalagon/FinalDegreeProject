'''
ALBA MALAGON MARQUEZ
Final Degree Project: 'STRATEGY FOR ACCURATE CNV DETECTION AND ML ALGORITHM FOR CLASSIFYING NGS VARIANTS'

Script similar to 'CNV_predictions.py' but with small modifications to analyse CNVs coming from arrays, not ExomeDepth. 

USAGE: python3 CNV_validation_method.py
'''

from CNV_libraries_variables import *


starttime = timeit.default_timer()


#    D E F I N I N G    V A R I A B L E S    ----------------------------------------------------------------------

thom, thom_, thet, thet_, fhom, fhom_, fhet, fhet_, incompdup, noinfo, dupnoinfo, compdup, sdups = 0,0,0,0,0,0,0,0,0,0,0,0,0
TP, FP = 0, 0
alldataGF = GFval  # all Genome Finals (SNVs and INDELs)
alldataED = EDvalcob # all Exome Depth files (CNVs) already covered


#    R E A D I N G    F I L E S   ----------------------------------------------------------------------------------

repeatmasker = pd.read_csv(RM[0].rsplit('/')[-1], sep="\t") # Repeatmasker 
pseudogenes = pd.read_csv(pseudogenes[0], sep='\t') # Pseudogenes 
segdups = pd.read_csv(SD[0], sep='\t') # Segmental duplications
truePdata = pd.read_csv(TPdata, sep="\t")
falsePdata = pd.read_csv(FPdata, sep="\t")


#    S T R A T E G Y    T O   P R E D I C T    C N V S   -----------------------------------------------------------


# creating two dictionaries to facilitate memorization
dict_GF={}
dict_ED={}

#adding a new column to the datasets
truePdata['info']='.'
falsePdata['info']='.'

#creating new datasets, which will contain the information the method predicts
truePdataM = pd.DataFrame(columns = truePdata.columns)
falsePdataM = pd.DataFrame(columns = falsePdata.columns)
noinfodataM = pd.DataFrame(columns = falsePdata.columns)
superdups_pseudoM = pd.DataFrame(columns = falsePdata.columns)
allcnvsM = pd.DataFrame(columns = truePdata.columns)
allcnvsM['num_variants']='.'
allcnvsM['prob_random']='.'

# storing the GENOME FINAL data in the dictionary in this way:    
# (key) qG20085488  -->  (value) qG20085488.genome.FINAL.txt
for f in alldataGF:
    identifier_GF = f.rsplit('/')[-1].rsplit('.')[0] 
    dict_GF[identifier_GF]=f


io=0 #a counter
for file_ED in alldataED:
    io+=1
    # storing the EXOME DEPTH data in the dictionary in this way:    
    # (key) qG20085488  -->  (value) qG20085488.ExomeDepthV3.PerSexe.female.2021_01_19.txt
    identifier_ED = file_ED.rsplit('/')[-1].rsplit('_')[0] 
    file_GF=dict_GF[identifier_ED]

    color_print('Analysing data {}     ({}/{}) ------------------------------------------------------------------'.format(identifier_ED, io, len(alldataED)), color='yellow')

    # reading the ExomeDepth and GenomeFinal files as dataframes using pandas
    exomedepth = pd.read_csv(file_ED, sep="\t")
    genomefinal = pd.read_csv(file_GF, sep="\t",dtype={"Chr":  "str"},low_memory=False)
    

    # getting the gender by looking at the chromosomes
    chry = exomedepth[(exomedepth['chr'].astype(str)=='Y')]
    if len(chry) > 0:
        gender_ED = 'male'
    else:
        gender_ED = 'female'
    

    # analysing each CNV individually
    for index, cnv in exomedepth.iterrows():
        chro = cnv['chr']
        start = cnv['start']
        end = cnv['end']
        typee = cnv['type']
        allchr='chr'+str(cnv['chr'])
        PARs,outsidePARs=False,False
        output1,output2=False,False
        repeatsM,SuperDups,partialSuperDup=False,False,False
        TP,FP,NI,SD_PS=False,False,False,False
        stillcompatible=False

        # taking only the variants falling within the CNV
        affected=genomefinal[(start<=genomefinal["Start"]) & (end>=genomefinal["End"]) & (chro==genomefinal["Chr"])]
        affected.reset_index(inplace=True)
        num_allvariants_affected = len(affected)
        allvariants=pd.DataFrame(affected)
        
        #SEGMENTAL DUPLICATIONS:  filtering CNVs overlapping with segmental duplications    
        sd1 = segdups[(start>=segdups["chromStart"]) & (end<=segdups["chromEnd"]) & (allchr==segdups["chrom"])]
        sd2 = segdups[(start<segdups["chromStart"]) & (end>segdups["chromEnd"]) & (allchr==segdups["chrom"])] #partial superdup
        sd3 = segdups[((start<segdups["chromStart"]) & (end>=segdups["chromStart"]) & (end<=segdups["chromEnd"]) & (allchr==segdups["chrom"]))] #partial superdup
        sd4 = segdups[((start<=segdups["chromEnd"]) & (start>=segdups["chromStart"]) & (end>segdups["chromEnd"]) & (allchr==segdups["chrom"]))] #partial superdup


        # these are the CNVs that fall completely into segmental duplications regions
        # since we cannot analyse them, we classify them as 'no info: whole cnv in segdup'
        if (len(sd1) > 0):

            SuperDups=True
            sdups+=1
            num='0'
            info = 'NO INFO: WHOLE CNV IN SEGDUP'
            c='cyan'
            output2=True

            sd2 = []
            sd3 = []
            sd4 = []


        # these are the CNVs that overlap partially with segmental duplications regions
        # we will analyze them, but eliminating previously the affected variants that overlap with the superdup

        if (len(sd2) > 0): 
            partialSuperDup=True
            if len(affected) > 0:
                for xs2, ys2 in affected.iterrows():
                    case2=sd2[(ys2['Start']>=sd2["chromStart"]) & (ys2['End']<=sd2["chromEnd"]) & (allchr==sd2["chrom"])] 
                    if len(case2) > 0:
                        affected.drop(xs2, axis=0, inplace=True)

        if len(sd3) > 0:
            partialSuperDup=True
            if len(affected) > 0:
                for xs3, ys3 in affected.iterrows():
                    case3=sd3[(ys3['Start']>=sd3["chromStart"]) & (ys3['End']<=sd3["chromEnd"]) & (allchr==sd3["chrom"])]
                    if len(case3) > 0:
                        affected.drop(xs3, axis=0, inplace=True)

        if (len(sd4) > 0):
            partialSuperDup=True
            if len(affected) > 0:
                for xs4, ys4 in affected.iterrows():
                    case4=sd4[(ys4['Start']>=sd4["chromStart"]) & (ys4['End']<=sd4["chromEnd"]) & (allchr==sd4["chrom"])]
                    if len(case4) > 0:
                        affected.drop(xs4, axis=0, inplace=True)


        #if the CNV does not fall completely within the segmental duplication, we can continue its analysis
        if SuperDups==False:       

            #REPEATMASKER:  filtering affected variants falling within repetitive regions (appearing in repeatmasker) 
            # also filtering the variants that do not have PASS in the Filter   
            if len(affected) > 0:
                filterNOTpass_index = affected[(affected['Filter'] != 'PASS')]
                affected = affected.drop(filterNOTpass_index.index, axis=0) 
                for xr, yr in affected.iterrows():
                    repeats=repeatmasker[(yr['Start']>=repeatmasker["genoStart"]) & (yr['End']<=repeatmasker["genoEnd"]) & (('chr'+str(yr['Chr']))==repeatmasker["genoName"])]
                    if len(repeats) > 0:
                        repeatsM=True
                        affected.drop(xr, axis=0, inplace=True)
                    
            #storing variables
            num_hom_a=len(affected[(affected['State']=='hom')])
            het_a=affected[(affected['State']=='het')]
            num_het_a=len(het_a)
            #taking the variants that have allele frequencies around 50
            freqAB=affected[(affected['AB']<0.65) | (affected['AB']>0.35)]    #50-50

            #PSEUDOAUTOSOMAL REGIONS
            #chrY:10001-2649520 and chrY:59034050-59363566
            #chrX:60001-2699520 and chrX:154931044-155260560
            startpar1 = 60001
            endpar1 = 2699520
            startpar2 = 154931044
            endpar2 = 155260560
            if allchr == 'chrX' and gender_ED=='male':
                PARs=True
                if (cnv['start']>endpar1) and cnv['end']<startpar2:
                    infoPARs= 'Considering a CNV of chrX of a male... and it falls OUTSIDE the PARs. So do not expect variants.'
                    outsidePARs=True
                else:
                    infoPARs= 'Considering a CNV of chrX of a male... but it falls INSIDE the PARs.'



            #considering the number of reads to know if we are in front of a deletion or a duplication

            if ((cnv['reads.observed']/cnv['reads.expected']) < 1): #deletion case
                typee='deletion'
                if ((cnv['reads.observed']/cnv['reads.expected']) <= 0.15) or (outsidePARs==True): #homozygous deletion case (or outside PARs)
                    if len(affected) >0: #homozygous deletion with variants
                        # how random is to find those variants there:
                        prob_var_random = 1
                        for x, y in affected.iterrows():
                            p = 0.02    # probability of success on a single trial (error)
                            n = y['DP']  # number of trials (depth)
                            k = (1-y['AB'])*y['DP'] # number of successes (number of alternatives)
                            var_randomly = binom.pmf(k, n, p)
                            prob_var_random *= var_randomly
                        if prob_var_random > min_prob_freq_random: #cannot trust the variants
                            num ='14'
                            info = 'COMPATIBLE HOMOZYGOUS DELETION'
                            thom+=1
                            output1=True
                            TP=True
                            c='green'
                            prob=prob_var_random

                        else: #trusting the variants

                            if cnv['reads.expected'] >= 100: 
                                num='16'
                                info = 'INCOMPATIBLE HOMOZYGOUS DELETION (has variants)'
                                fhom+=1
                                output1=True
                                FP=True
                                c='red'
                                prob=prob_var_random
                            
                            else:
                                num='15'
                                info = 'INCOMPATIBLE HOMOZYGOUS DELETION... CHECK! (has variants)'
                                fhom_ += 1
                                output1=True
                                FP=True
                                c='red'
                                prob=prob_var_random

                    else: ##homozygous deletion without variants
                        if cnv['reads.expected'] >= 100: 
                            num='13'
                            info = 'COMPATIBLE HOMOZYGOUS DELETION'
                            thom+=1
                            output2=True
                            TP=True
                            c='green'

                        else:
                            num='12'
                            info = 'COMPATIBLE HOMOZYGOUS DELETION... CHECK!'
                            thom_+=1
                            output2=True
                            TP=True
                            c='green'

                        
                else: #heterozygous deletion
                    if len(affected)==0: #heterozygous deletion without variants. No variants means no info
                        num='17'
                        info = 'NO INFO DELETION'
                        noinfo+=1
                        output2=True
                        c='magenta'
                        NI=True

                    else:  #heterozygous deletion with variants                  
                        if num_het_a > 0: #at least there is an heterozygous variant within the CNV
                            # assuming an error of 2% (P)... what is the probability of observing K alternatives with N depth? binomial probability
                            prob_het_random = 1
                            for i, h in het_a.iterrows():
                                p = 0.02    # probability of success on a single trial (error)
                                n = h['DP']  # number of trials (depth)
                                #print(h)
                                k = (1-h['AB'])*h['DP'] # number of successes (number of alternatives)
                                het_randomly = binom.pmf(k, n, p)
                                prob_het_random *= het_randomly

                            if prob_het_random < min_prob_freq_random: #trusting the variants
                                num='11'
                                info = 'INCOMPATIBLE HETEROZYGOUS DELETION (has het)'
                                fhet+=1
                                output1=True
                                FP=True
                                c='red'
                                prob=prob_het_random
                            
                            else: #since we cannot trust the het variants, it is still compatible, so we continue the analysis
                                stillcompatible=True


                        # considering the case where all affected variants are HOM (or that the het variants detected are not trusted)
                        elif (num_het_a<=0) or (stillcompatible==True):
                            # if there are not enough variants affecting the region (< 3), it would make no sense to look for controls
                            # the CNVs that does not satisfy this requirement will be filtered as 'no info'
                            if len(affected) < min_variants_affected:
                                num='7'
                                info = 'NO INFO DELETION'
                                noinfo+=1
                                output2=True
                                c='magenta'
                                NI=True


                            # considering the presence of controls
                            else:
                                counting_variants_control = 0
                                num_hom_c = 0
                                num_het_c = 0
                                ab = []
                                
                                # generating a dataframe with all ED files, excluding the reference test
                                all_ED = pd.DataFrame(columns = exomedepth.columns)
                                c_ed = [ i for i in alldataED if i!=file_ED ]   
                                for control_ed in c_ed:
                                    ED=pd.read_csv(control_ed, sep="\t")
                                    all_ED = all_ED.append(ED)

                                # generating a dataframe with all GF files, excluding the reference test
                                controls = pd.DataFrame(columns = genomefinal.columns)
                                c_gf = [ i for i in alldataGF if i!=file_GF]
                                for control_gf in c_gf:
                                    GF=pd.read_csv(control_gf, sep="\t",dtype={"Chr":  "str"},low_memory=False)
                                    chroo = cnv['chr']
                                    startt = cnv['start']
                                    endd = cnv['end']
                                    # the possible control variants are the ones within the cnv region of the reference
                                    possible_c=GF[(startt<=GF["Start"]) & (endd>=GF["End"]) & (chroo==GF["Chr"])]
                                    # for each of these possible control variants, we need to ensure that they do not belong to any other CNV of any other ED file
                                    for index2, candidate in possible_c.iterrows():
                                        startt2=candidate["Start"]
                                        endd2=candidate["End"]
                                        chroo2=candidate["Chr"]
                                        if len(all_ED[(startt2>=all_ED["start"]) & (endd2<=all_ED["end"]) & (str(chroo2)>=str(all_ED["chr"]))])==0:
                                            controls=controls.append(candidate)

                                counting_variants_control=len(controls)
                                num_hom_c+=len(controls[(controls["State"]=="hom")])
                                num_het_c+=len(controls[(controls["State"]=="het")])
                                ab.append(controls['AB'])
                            

                                #if there are not enough control variants in the region (< 3), the region cannot be compared with the null distribution.
                                #the CNV that does not satisfy this requirement will be filtered as 'still compatible heterozygous deletion'

                                if counting_variants_control < min_variants_control:
                                    num='8'
                                    info = 'still COMPATIBLE HETEROZYGOUS DELETION... CHECK'
                                    thet_ +=1
                                    output2=True
                                    TP=True
                                    c='green'

                                else:
                                    num='9'
                                    info = 'COMPATIBLE HETEROZYGOUS DELETION (all hom)'
                                    thet+=1
                                    output1=True
                                    TP=True
                                    c='green'
                                    # how strange is that all variant are HOM when considering a heterozygous:
                                    binomial_prob = binom.pmf(k=num_het_a, n=len(affected), p=num_het_c/counting_variants_control)
                                    prob=binomial_prob

                                    print()
                                    plt.hist(ab, density=True, color='cyan')
                                    plt.xlabel('allele balance')
                                    plt.ylabel('frequency') 
                                    plt.title('{} cnv{:n}'.format(identifier_ED, index+1), color='cyan') 
                                    plt.suptitle('control AB distribution', color='cyan')
                                    plt.savefig('plotAB2_{}_cnv{:n}.png'.format(identifier_ED, index+1))
                                    plt.close() 


            else: #duplication case
                typee='duplication'
                if (cnv['reads.observed']/cnv['reads.expected']>1.2) and (cnv['reads.observed']/cnv['reads.expected']<1.7): #three copies
                    if (num_het_a>0): #it has at least one heterozygous variants
                        if len(freqAB) > 0: # variants with 50-50
                            #the null hypothesis is that I am trisomic… so the freq is 30%… what is the proof that I can see 50-50?
                            # how random is to find those variables there:
                            prob_freq_random = 1
                            for x, y in freqAB.iterrows():
                                p = 0.30    # probability of success on a single trial 
                                n = y['DP']  # number of trials (depth)
                                k = (1-y['AB'])*y['DP'] # number of successes (number of alternatives)
                                freq_randomly = binom.pmf(k, n, p)
                                prob_freq_random *= freq_randomly
                            if prob_freq_random < min_prob_freq_random: #trusting the variants
                                num='6'
                                info = 'INCOMPATIBLE DUPLICATION (3 copies and AB 50-50)'
                                incompdup+=1
                                output1=True
                                FP = True
                                c='yellow'
                                prob=prob_freq_random

                            else: #cannot trust the variants
                                num='5'
                                info = 'COMPATIBLE DUPLICATION'
                                compdup+=1
                                output1=True
                                TP = True
                                c='cyan'
                                prob=prob_freq_random

                    
                        else: #3 copies and AB (<0.4 and >0.6)
                            # how random is to find those variables there:
                            prob_freq_random = 1
                            for x, y in freqAB.iterrows():
                                p = 0.50    # probability of success on a single trial 
                                n = y['DP']  # number of trials (depth)
                                k = (1-y['AB'])*y['DP'] # number of successes (number of alternatives)
                                freq_randomly = binom.pmf(k, n, p)
                                prob_freq_random *= freq_randomly
                            if prob_freq_random < min_prob_freq_random: #trusting the variants
                                num='4'
                                info = 'COMPATIBLE DUPLICATION (3 copies and AB (<0.4 and >0.6))'
                                compdup+=1
                                output1=True
                                TP = True
                                c='cyan'
                                prob=prob_freq_random

                            else: #cannot trust the variants
                                num='3'
                                info = 'INCOMPATIBLE DUPLICATION'
                                incompdup+=1
                                output1=True
                                FP = True
                                c='yellow'
                                prob=prob_freq_random

                    else:  # no heterozygous variants to analyse
                        num='2'
                        info = 'NO INFO DUPLICATION (3 copies but no variants)'
                        dupnoinfo+=1
                        output2=True
                        NI = True
                        c='magenta'

                else: #duplication with more than 3 copies
                    num='1'
                    info = 'NO INFO DUPLICATION (+3 copies)'
                    dupnoinfo+=1
                    output2 = True
                    NI = True
                    c='magenta'



        # printing the output
        print()
        if partialSuperDup:
            info += ' partial superdup'

        print()
        color_print('{} --------------------------------------- {} ---------------------------------------'.format(num,info), color=c)
        print()

        if output1: #it contains the probability of the variants appearing at random
            print('EXOME DEPTH INFORMATION:')
            color_print(('\tType: {}\n\tID: {}\n\tCNV: {}\n\tChr: {}\n\tStart: {}\n\tEnd: {}\n\tobserved reads: {}\n\texpected reads: {}\n\tratio {}\n\tall variants: {}\n\tsignificant variants: {}\n\tBF: {}\n\tProbability at random: {}\n\tGender: {}'.format(typee,identifier_ED,index+1, chro, start, end, cnv['reads.observed'], cnv['reads.expected'], cnv['reads.observed']/cnv['reads.expected'], num_allvariants_affected, len(affected), cnv['BF'], prob, gender_ED)), color='white')
        if output2: #no probability score
            print('EXOME DEPTH INFORMATION:')
            color_print(('\tType: {}\n\tID: {}\n\tCNV: {}\n\tChr: {}\n\tStart: {}\n\tEnd: {}\n\tobserved reads: {}\n\texpected reads: {}\n\tratio {}\n\tall variants: {}\n\tsignificant variants: {}\n\tBF: {}\n\tGender: {}'.format(typee,identifier_ED, index+1, chro, start, end, cnv['reads.observed'], cnv['reads.expected'], cnv['reads.observed']/cnv['reads.expected'], num_allvariants_affected, len(affected), cnv['BF'], gender_ED)),color='white')
        if len(affected) > 0 and SuperDups == False:
            print()
            print('SIGNIFICANT VARIANTS WITHIN THE CNV:') #printing the significant variants!
            print(affected[['Start','End','Chr','Ref','Alt','State','AB','DP','Filter','Func.refGene','Gene.refGene','genomicSuperDups']] )
            print()
            # PSEUDOGENES
            pseudo=pseudogenes[(pseudogenes['cdsStart']>=cnv['start']) & (pseudogenes['cdsEnd']<=cnv['end']) & (pseudogenes['#chrom']==allchr)]
            pseudo.reset_index(inplace=True)
            if len(pseudo) > 0:    
                SD_PS=True
                info += ' - PSEUDOGENE'
                print('Be careful, the type of gene we are considering is a pseudogene.\n')
        
        if PARs: #info PARs
            print('\n{}\n'.format(infoPARs))

        if repeatsM: #info RepeatMasker
            print()
            color_print('{}'.format('The following variants were falling within a repetitive region, so they have been excluded'),color='white')
            color_print('\nVariant indentified      -->     start:{}   end:{}   chr:{}  state:{}    AB:{}   DP:{}   Ref:{}  Alt:{}  Gene:{}'.format(yr['Start'], yr['End'], yr['Chr'],yr['State'],yr['AB'],yr['DP'],yr['Ref'],yr['Alt'],yr['Gene.refGene']), color='white')
            print('\nRepeatMasker information:')
            print(repeats)

        
        if SuperDups: #info Segmental Duplications
            color_print('{}'.format('\nWhole CNV falling within a SuperDup!!!'),color='cyan')
            SD_PS=True
            if len(sd1) > 0:
                print('CASE 1: complete CNV falling within a segmental duplication')
                print(sd1)

            if num_allvariants_affected > 0:
                print()
                print('ALL VARIANTS WITHIN THE CNV:') #printing all variants: significant or not
                print(affected[['Start','End','Chr','Ref','Alt','State','AB','DP','Filter','Func.refGene','Gene.refGene','genomicSuperDups']] )
                print()

        if partialSuperDup: #info partial Segmental Duplications
            color_print('{}'.format('\nCNV falling partially with SuperDups!!!'),color='yellow')
            if len(sd2) > 0:
                print('CASE 2')
                print(sd2)

            if len(sd3) > 0:
                print('CASE 3')
                print(sd3)
            
            if len(sd4) > 0:
                print('CASE 4')
                print(sd4)

            if num_allvariants_affected > 0:
                print()
                print('ALL VARIANTS WITHIN THE CNV:') #printing all variants: significant or not
                print(allvariants[['Start','End','Chr','Ref','Alt','State','AB','DP','Filter','Func.refGene','Gene.refGene','genomicSuperDups']] )
                print()

        if SD_PS:
            superdups_pseudoM=superdups_pseudoM.append(cnv)
            superdups_pseudoM['info'].iloc[-1] =info
        
        else:   #do not contain a superdup nor a pseudogene
            if TP:
                truePdataM = truePdataM.append(cnv)
                truePdataM['info'].iloc[-1] =info

            if FP:
                falsePdataM = falsePdataM.append(cnv)
                falsePdataM['info'].iloc[-1] =info

            if NI:
                noinfodataM = noinfodataM.append(cnv)
                noinfodataM['info'].iloc[-1] =info


        allcnvsM=allcnvsM.append(cnv)
        allcnvsM['info'].iloc[-1]=info
        allcnvsM['num_variants'].iloc[-1]=len(affected)
        if output1:
            allcnvsM['prob_random'].iloc[-1]=prob

#removing features 
truePdataM.drop(['OMIM', 'OMIMFeb2017', 'NqExDels',
       'NqExDups', 'qExDels', 'qExDups', 'Unnamed: 0'], axis=1, inplace=True)
falsePdataM.drop(['OMIM', 'OMIMFeb2017', 'NqExDels',
       'NqExDups', 'qExDels', 'qExDups', 'Unnamed: 0'], axis=1, inplace=True)
noinfodataM.drop(['OMIM', 'OMIMFeb2017', 'NqExDels',
       'NqExDups', 'qExDels', 'qExDups', 'Unnamed: 0'], axis=1, inplace=True)
superdups_pseudoM.drop(['OMIM', 'OMIMFeb2017', 'NqExDels',
       'NqExDups', 'qExDels', 'qExDups', 'Unnamed: 0'], axis=1, inplace=True)

#printing general CNVs' counts
color_print('COMPATIBLE:     {}'.format(len(truePdataM)), color='green')
color_print('NO COMPATIBLE:     {}'.format(len(falsePdataM)), color='red')
color_print('NO INFO:     {}'.format(len(noinfodataM)), color='yellow')
color_print('SuperDups and pseudogenes:     {}'.format(len(superdups_pseudoM)), color='magenta')

#storing datasets
truePdataM.to_csv('{}/truePdataMETHOD.csv'.format(os.getcwd()), header=True, sep="\t")
falsePdataM.to_csv('{}/falsePdataMETHOD.csv'.format(os.getcwd()), header=True, sep="\t")
noinfodataM.to_csv('{}/NOINFOdataMETHOD.csv'.format(os.getcwd()), header=True, sep="\t")
superdups_pseudoM.to_csv('{}/superdups_pseudoMETHOD.csv'.format(os.getcwd()), header=True, sep="\t")
allcnvsM.to_csv('{}/allcnvsMETHOD.csv'.format(os.getcwd()), header=True, sep="\t")

#computing compatible (TP) and incompatible (FP) CNVs
TP = thom + thom_ + thet + thet_
FP = fhom + fhom_ + fhet + fhet_

                
# NORMAL PIE CHART
labels=['INCOMPATIBLE HOMOZYGOUS DELETION','COMPATIBLE HOMOZYGOUS DELETION', 'COMPATIBLE HETEROZYGOUS DELETION','COMPATIBLE HOMOZYGOUS DELETION...CHECK', 'INCOMPATIBLE HOMOZYGOUS DELETION...CHECK', 'INCOMPATIBLE HETEROZYGOUS DELETION', 'COMPATIBLE HETEROZYGOUS DELETION...CHECK','INCOMPATIBLE DUPLICATION', 'INCOMPATIBLE HETEROZYGOUS DELETION...CHECK','NO INFO', 'DUP NO INFO', 'COMPATIBLE DUPLICATION', 'SUPERDUPS']
sizes=[fhom,thom,thet,thom_,fhom_,fhet,thet_, incompdup, fhet_,noinfo, dupnoinfo, compdup, sdups]
colors = ['magenta', 'yellow', 'purple', 'green', 'red', 'orange', 'cyan', 'pink', 'black', 'white', 'blue', 'grey','cyan']
totalcnvs = sum(sizes)
percentages=[]
for each in sizes:
    percentages.append(str(round(((each/totalcnvs)*100),2))+'%')
plt.pie(sizes, colors=colors, labels=percentages, shadow=True, startangle=90)
plt.legend( loc = 'lower right', labels=labels, title='SIMPLE CNVs PIE CHART')
plt.axis('equal')
plt.tight_layout()
plt.savefig('PIE_validation.png')
plt.close()



# PLOT COMPATIBILITY
labels = ['compatible deletion', 'incompatible deletion', 'false duplication', 'no info', 'dup no info', 'true dup', 'superdups']
sizes = [TP, FP, incompdup, noinfo, dupnoinfo, compdup,sdups]
colors = ['magenta', 'cyan', 'yellow', 'white', 'purple', 'orange','red']
plt.pie(sizes, colors=colors, autopct='%1.1f%%', shadow=True, startangle=90)
plt.legend( loc = 'lower left', labels=labels, title='CNVs COMPATIBILITY PIE CHART')
plt.axis('equal')
plt.tight_layout()
plt.savefig('PIEcompatibility_validation.png')
plt.close()



# NESTED PIE CHART
group_names=['COMPATIBLE HOMOZYGOUS DELETION','INCOMPATIBLE HOMOZYGOUS DELETION', 'INCOMPATIBLE HETEROZYGOUS DELETION', 'COMPATIBLE HETEROZYGOUS DELETION', 'DUPLICATION', 'NO INFO']
group_size=[thom+thom_,fhom+fhom_,fhet+fhet_,thet+thet_, incompdup+dupnoinfo+compdup, noinfo+sdups]
subgroup_names=['COMPATIBLE HOMOZYGOUS DELETION', 'COMPATIBLE HOMOZYGOUS DELETION...CHECK', 'INCOMPATIBLE HOMOZYGOUS DELETION', 'INCOMPATIBLE HOMOZYGOUS DELETION...CHECK', 'INCOMPATIBLE HETEROZYGOUS DELETION', 'INCOMPATIBLE HETEROZYGOUS DELETION...CHECK', 'COMPATIBLE HETEROZYGOUS DELETION', 'COMPATIBLE HETEROZYGOUS DELETION...CHECK', 'DUPLICATION', 'DUP NO INFO', 'DUP GOOD', 'NO INFO', 'SUPERDUPS']
subgroup_size=[thom,thom_,fhom,fhom_,fhet, fhet_,thet,thet_, incompdup, dupnoinfo, compdup, noinfo, sdups]
totalcnvs = sum(group_size)
percentages=[]
for each in group_size:
    percentages.append(str(round(((each/totalcnvs)*100),2))+'%')

# Create colors
one, two, three, four, five, six=[plt.cm.Blues, plt.cm.Reds, plt.cm.Greens, plt.cm.Purples, plt.cm.Oranges, plt.cm.Greys]

# First Ring (outside)
fig, ax = plt.subplots(figsize =(10, 7))
ax.axis('equal')
mypie, _ = ax.pie(group_size, radius=1.1, labels=percentages, colors= 
[one(0.6), two(0.6), three(0.6), four(0.6), five(0.6), six(0.1)] )
plt.setp( mypie, width=0.5, edgecolor='white')

# Second Ring (Inside)
mypie2, _ = ax.pie(subgroup_size, radius=1.1-0.5, 
colors=[one(0.5), one(0.4), 
two(0.3), two(0.5), three(0.4), three(0.6), four(0.5), four(0.3), five(0.6), five(0.3), five(0.4),six(0.1),six(0.4)])
plt.setp( mypie2, width=0.6, edgecolor='white')
ax.set_title("CNVs NESTED PIE CHART")
ax.legend(group_names, loc='lower left')

plt.savefig('NestedPIE_validation.png')
plt.close()






stoptime = timeit.default_timer()
print('\n\nTime:',stoptime - starttime)




