import subprocess  # To run commands in the terminal
import os  # To navigate through system directories
import click  # Library to create command-line commands
@click.command()
@click.option('--taxon-id', prompt='Enter the taxon ID (Streptomyces: 1883)',
              help='The taxon ID to be used in the command.')


def run_datasets_summary(taxon_id): 
  
    # Build the command
    command = ['./datasets', 'download', 'genome','gbff', 'taxon', str(taxon_id), '--annotated' '--seq-report']
    #Options to the command '--include','protein,genome,rna,cds,gff3,gtf,gbff,seq-report']
   
    # Print the command for debugging purposes
    print(f"Running command: {' '.join(command)}")

    try:
         #Run the command
        subprocess.run(command, check=True)
        #After the command is run, move or store the result in the download folder
        print(f"Command 'datasets download' executed successfully for taxon {taxon_id}.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing the command: {e}")
    except FileNotFoundError:
       print('Error: The "./datasets" file was not found. Please check if the path is correct.')

if __name__ == '__main__':
    run_datasets_summary()