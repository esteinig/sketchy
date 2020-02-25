import click
import psutil
import time
import logging
import pandas

from sketchy.utils import PoreLogger


@click.command()
@click.option(
    '--interval', '-i', default=0.1, type=float, required=False,
    help='Interval in seconds to check on resource usage'
)
@click.option(
    '--terminate', '-t', is_flag=True, required=False,
    help='Terminate monitoring if no more processes are active'
)
@click.option(
    '--early', '-e', default=None, type=float, required=False,
    help='Terminate early after this amount of seconds'
)
@click.option(
    '--prefix', '-p', default='sketchy.benchmark', type=str,
    help='Prefix for data output file when completed or stopped',
)
@click.option(
    '--log', '-l', is_flag=True,
    help='Output to log file instead of terminal'
)
def monitor(interval, terminate, early, prefix, log):

    """ Monitor benchmarks during a Sketchy execution (Mash, Sketchy) """

    # A bit hacky because Sketchy calls Mash in parallel, and there is
    # currently no return of specific PID from Mash

    # Essentially monitors ANY command that contains: mash, sketchy-rs
    # If multiple commands are running, resource usage is computed across all

    if interval < 0.1:
        raise ValueError('Interval (--interval, -i) must be >= 0.1')

    start_time = time.time()

    if log:
        logfile = f'{prefix}.log'
    else:
        logfile = None

    outfile = f'{prefix}.tsv'

    pore = PoreLogger(level=logging.INFO, file=logfile)

    data = []

    try:
        pore.logger.info(f'CPU\tMEM\tMash\t\tSketchy\t')
        pore.logger.info(f'---\t---\t----\t\t-------\t')
        while True:
            run_time = time.time() - start_time
            if early and run_time > early:
                pore.logger.info(f'Early termination after {terminate} seconds.')
                exit(0)

            sketchy_pids, mash_pids = check_pids()

            # Early termination if no processes
            if not sketchy_pids and not mash_pids:
                if terminate:
                    pore.logger.info('No more processes found. Exiting.')
                    summarize_data(data, outfile, pore.logger)
                    exit(0)
                else:
                    pass
            else:
                try:
                    cpu1, rss1 = get_resource_use(mash_pids)
                    cpu2, rss2 = get_resource_use(sketchy_pids)

                    cpu_total = cpu1 + cpu2
                    rss_total = rss1 + rss2

                    pore.logger.info(
                        f'{round(cpu_total, 1)}\t{round(rss_total, 1)}\t'
                        f'{round(cpu1, 1)}\t{round(rss1, 1)}\t'
                        f'{cpu2}\t{round(rss2, 1)}'
                    )

                    data.append(
                        [run_time, cpu_total, rss_total, cpu1, rss1, cpu2, rss2]
                    )
                except psutil.NoSuchProcess:
                    if terminate:
                        pore.logger.info('No more processes found. Exiting.')
                        summarize_data(data, outfile, pore.logger)
                        exit(0)
                    time.sleep(interval)

            time.sleep(interval)

    except KeyboardInterrupt:
        summarize_data(data, outfile, pore.logger)
        exit(0)


def summarize_data(data, fout, logger):

    df = pandas.DataFrame(
        data, columns=[
            'time', 'cpu', 'mem', 'cpu_mash', 'mem_mash',
            'cpu_sketchy', 'mem_sketchy'
        ]
    )

    df = df[df.cpu != 0]  # sometimes happens to catch Rust library on init

    mean, std, peak = df.mean(), df.std(), df.max()

    mean_cpu = round(mean["cpu"], 2)
    mean_cpu_std = round(std["cpu"], 2)
    peak_cpu = round(peak["cpu"], 2)

    mean_mem = round(mean["mem"], 2)
    mean_mem_std = round(std["mem"], 2)
    peak_mem = round(peak["mem"], 2)

    logger.info(
        f'CPU: {mean_cpu}% +- {mean_cpu_std}% (max: {peak_cpu}%)'
    )
    logger.info(
        f'Memory: {mean_mem} MB +- {mean_mem_std} MB (max: {peak_mem} MB)'
    )

    df.to_csv(fout, sep='\t', index=False)

    print(
        f'{mean_cpu}\t{mean_cpu_std}\t{peak_cpu}\t'
        f'{mean_mem}\t{mean_mem_std}\t{peak_mem}'
    )


def get_resource_use(pids: list) -> tuple:

    total_cpu, total_rss = 0, 0
    for pid in pids:
        p = psutil.Process(pid)
        total_cpu += p.cpu_percent(0.1)
        total_rss += p.memory_info().rss / float(1 << 20)

    return total_cpu, total_rss


def check_pids():

    sketchy_pids = []
    mash_pids = []

    for process in psutil.process_iter():
        cmdline = process.cmdline()
        if 'sketchy-rs' in cmdline:
            sketchy_pids.append(process.pid)
        if 'mash' in cmdline:
            mash_pids.append(process.pid)

    return sketchy_pids, mash_pids
