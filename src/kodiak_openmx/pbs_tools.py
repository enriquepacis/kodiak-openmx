from datetime import datetime
import subprocess as sp
import os, re

def start_status_file( fname, job_id ):
    """Initializes a PBS status file, recoding the
    PBS job ID and the start time.
    """

    now = datetime.now()
    time_str = now.strftime("%b %d %H:%M:%S %Y")

    with open(fname, 'a') as statf:
        statf.write(f'Job queued ({job_id}). {time_str}\n')

def get_job_id( status_file ):

    with open( status_file ) as f:
        text = f.read()

    match = re.search(r"\((\d+)\)", text)

    if match:
        job_id = match.group(1)
        print(f"Job ID: {job_id}")
    else:
        job_id = None
        print("No job ID found.")

    return job_id

def job_done(outfile):
    """
       This queries the output file to determine whether
       the job is done. The output file must exist and
       contain the key string 'no' for this method
       to return a Python True value.
    """

    jobdone = False

    # assert os.path.isfile( outfile ), f'Output file ({outfile}) does not exist.'

    if os.path.isfile( outfile ):
        with open( outfile, 'r') as f:
            outf_lines = f.readlines()

        for line in outf_lines:
            jobdone = jobdone or 'normally finished' in line

    return jobdone

def get_job_status( status_file, output_file ):
    """
       Obtains the status of a job given a status file name
       and an output file name.

       This is designed for non-pw.x calculations, such as
       dos.x or pw.x

       SYNTAX
       ======
       job_status = get_job_status( status_file, output_file)

       if job_status == 'done':
           <do some post-processing>
    """

    status = None
    JobDone = False

    hasStatus = os.path.isfile(status_file)
    hasOutput = os.path.isfile(output_file)

    # assert FoundStatus or FoundLocalStatus, 'Status file not found.'
    # print(os.getcwd())
    # print(f'status file: {status_file}\noutput file: {output_file}')
    # print('Status file exists: {0}'.format( os.path.isfile( status_file ) ))
    # print('Output file exists: {0}'.format( os.path.isfile( output_file ) ))
    if hasStatus:

        # Extract lines of status file
        with open(status_file, 'r') as sf:
            status_lines = sf.readlines()

        # for line in status_lines:
        #    print(line)

        # Analyze lines of status file
        JobQueued = False
        JobStarted = False
        JobEnded = False
        JobID = None
        for line in status_lines:
            if 'Job queued' in line:
                JobQueued = True
                # print('Queue data: {0}'.format(line[:-1]))
                JobID = re.findall(r'\((.*?)\)', line)[0]
                # print(JobID)
                # print(f' Found JobID: {JobID}')
            if 'Commenced job' in line:
                JobStarted = True
            if 'Job ended' in line:
                JobEnded = True

        # check for JOB DONE status
        # print(f'JobQueued = {JobQueued}')
        # print(f'JobStarted = {JobStarted}')
        # print(f'JobEnded = {JobEnded}')

        if JobEnded:
            if job_done( output_file ):
                status = 'done'
            else:
                status = 'aborted'

        else: # job isn't done. It could be aborted, running, or queued

            if JobID is not None:
                # print(f'Querying qstat for job {JobID}.')
                status = PBS_job_status( 'scf_status.txt' )

            else:
                if JobQueued and not JobStarted:
                    status = 'queued'

                elif JobQueued and JobStarted:
                    status = 'running'
                    # check for aborted status

    elif hasOutput and not hasStatus:
        # Older calculations may not have a status file
        # But an output file is required
        if job_done( output_file ):
            status = 'done'
        else:       # Job is NOT DONE
            status = 'unknown'

    elif not (hasOutput or hasStatus):
        status = None

    return status

def PBS_job_status( status_file ):
    """
       Queries PBS and returns the PBS status of a job.

       If the job return value is None, the job
       is either done, aborted, or does not exists.
    """

    status = None

    jobID = get_job_id( status_file )

    qstat = sp.Popen('/bin/bash',
                     stdin=sp.PIPE,
                     stdout=sp.PIPE,
                     stderr=sp.PIPE,
                     universal_newlines=True)

    cmd_str = f'qstat -J {jobID}'
    out, err = qstat.communicate( cmd_str )

    out_lines = out.split('\n') # split the output by newlines

    # print('The output has {0} lines: "{1}"'.format(len(out_lines), out))

    """
       If the job jobID doesn't exist, out = '\n', so this counts
           as one line after the out.split('\n').

       If the job exists, then 'qstat -j {jobID}' returns a header
       line, a divider line, and a data line (plus carriage return)
       for a total of 4 lines.
    """
    if len(out_lines) > 1:
        """

        for line in out_lines:
            print(line)
        """
        status_line = out_lines[2] # data line
        status_char = list(filter(None, status_line.split()))[4]

        # print(f'status_char = {status_char}')

        status_dict = {'R': 'running', 'Q': 'queued', 'E': 'error',
                       'A': 'aborting'}

        status = status_dict[status_char]

    return status
