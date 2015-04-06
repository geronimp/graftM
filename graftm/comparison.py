        # If there are reverse reads
        if len(reads_list) == 2:
            Messenger().message('Found %s forward read(s) for %s' % (len(reads_list[0]), base))
            Messenger().message('Found %s reverse read(s) for %s' % (len(reads_list[1]), base))
            write_list = [x for x in reads_list[1] if x in set(reads_list[0])]
            if len(write_list) > 0:
                Messenger().message('%s read(s) overlapped, continuing with them' % (len(write_list)))
            else:
                Messenger().message('%s read(s) overlapped, cannot continue' % (len(write_list)))
                exit(0)
