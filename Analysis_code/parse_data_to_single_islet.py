import numpy 
import glob

data_dirr = 'DevData/'

dev_stages = [\
            'Gestation/'\
            ,'1-35weeks/'\
            ,'12-24months/'
            ,'28MandAfter/'\
                    ]

dev_labels = [\
        '0'\
        ,'1'\
        ,'2'\
        ,'3'\
            ]

in_dirr = 'all/'

all_islets_info_file = 'all_islets_info.csv'
mm = open(all_islets_info_file, 'w')
mm.write('Stage,Subject,Islet,fileprefix,Area,n_ad,n_b\n')

#count = 0

prev_area = 0

for idx, dev in enumerate(dev_stages):

    this_dirr = data_dirr + dev + in_dirr

    all_subjects = glob.glob(this_dirr+'*.tsv')

    for subject in all_subjects:

        line = subject.split('/')
        subject_name = line[-1].split('.')[0]

        this_stage = dev_labels[idx]

        islet_prefix = 'Stage'+ this_stage + '_' + subject_name 

        ff = open(subject, 'r')

        isletnum = -1
        
        write_lines = []
        n_ad = 0
        n_b = 0

        #print(count, end='\r')

        for line in ff:
            this_line = line.split('\t')
            this_islet_num = int(this_line[0])
            ty = int(this_line[-2])

            this_area = float(this_line[-1])

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
                           + ',' + str(prev_area)\
                           + ',' + str(n_ad)\
                           + ',' + str(n_b)\
                           + '\n'
                           )

                        #gg = open(data_dirr + '/all_islets/' + isletfile + '.tsv', 'w')

                        #for lline in write_lines:
                        #    gg.write(lline)

                        #gg.close()

                    n_b = 0
                    n_ad = 0
                    write_lines = []

                isletnum = this_islet_num
                prev_area = this_area

            write_lines.append(line)

        ff.close()



mm.close()



