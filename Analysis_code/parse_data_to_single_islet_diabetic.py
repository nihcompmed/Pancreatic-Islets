
import numpy 
import glob

data_dirr = 'DiabeticData/'

stages = [\
            'Control/'\
            ,'Diabetic/'\
                    ]

labels = [\
        'C'\
        ,'D'\
            ]

all_islets_info_file = 'all_islets_info_diabetic.csv'
mm = open(all_islets_info_file, 'w')
mm.write('Stage,Subject,Islet,fileprefix,n_ad,n_b\n')

#count = 0

dummy_area = -1

for idx, stage in enumerate(stages):

    this_dirr = data_dirr + stage

    all_subjects = glob.glob(this_dirr+'*.tsv')
    # gfp is beta
    # rfp is a
    # cy5 is d

    for subject in all_subjects:

        line = subject.split('/')
        subject_name = line[-1].split('.')[0]

        this_stage = labels[idx]

        islet_prefix = 'Stage'+ this_stage + '_' + subject_name 

        ff = open(subject, 'r')

        isletnum = -1
        
        write_lines = []
        n_ad = 0
        n_b = 0

        #print(count, end='\r')

        for line in ff:

            line = line.strip('\n')
            if line == '':
                continue

            this_line = line.split('\t')
            #print(subject, this_line)

            this_islet_num = int(this_line[0])

            if this_line[-1] == 'rfp':
                ty = 1
            elif this_line[-1] == 'gfp':
                ty = 2
            elif this_line[-1] == 'cy5':
                ty = 3
            else:
                print(this_line[-1])
                exit()

            #this_area = float(this_line[-1])

            if ty == 2:
                n_b += 1
            else:
                n_ad += 1

            if isletnum != this_islet_num:

                if isletnum != -1:
                    if (n_ad > 5) and (n_b > 5):

                        isletfile = islet_prefix + '_Islet' + str(isletnum)

                        #('Stage,Subject,Islet,fileprefix,Area,n_ad,n_b\n')
                        mm.write(this_stage\
                           + ',' + subject_name\
                           + ',' + str(isletnum)\
                           + ',' + isletfile\
                           + ',' + str(dummy_area)\
                           + ',' + str(n_ad)\
                           + ',' + str(n_b)\
                           + '\n'
                           )

                        gg = open(data_dirr + '/all_islets/' + isletfile + '.tsv', 'w')


                        for lline in write_lines:
                            gg.write(lline)

                        gg.close()

                    n_b = 0
                    n_ad = 0
                    write_lines = []

                isletnum = this_islet_num

            edit_line = line.strip('\n')
            edit_line = line.split('\t')
            new_line = edit_line[0]\
                    + '\t' + edit_line[1]\
                    + '\t' + edit_line[2]\
                    + '\t' + edit_line[3]\
                    + '\t' + str(ty)\
                    + '\t' + str(dummy_area)\
                    + '\n'
            write_lines.append(new_line)

        ff.close()



mm.close()



