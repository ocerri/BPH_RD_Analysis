import os, psutil

# To get process:
#     process = psutil.Process(os.getpid())
def print_memory_usage(process, tag=''):
    mem = process.memory_info()
    B2MB = 9.537e-7
    B2GB = 9.31e-10
    print '\n' + 40*'#' + '\n' + 20*'#'
    print tag
    print 10*'-' + '> Mem: {:.1f} GB (res: {:.1f} GB, virt: {:.1f} GB)     Data Mem: {:.1f} MB'.format(B2GB*(mem.rss+mem.vms), B2GB*mem.rss, B2GB*mem.vms, B2MB*mem.data)
    print '\n' + 20*'#' + '\n' + 40*'#' + '\n'
