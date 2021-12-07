
import gnomadIC
import argparse
import hail as hl

def main(args):
    '''Controls whether to setup in test mode or not, and generates a run ID if not in test mode'''
    # Initialise Hail, setting output 
    hl.init(
        log='hail_logs/log.txt', 
        quiet=args.quiet
        )

    # Setup paths
    if args.test:
        run_ID = 'test'
        print('Running in test mode: Relax, sit back and enjoy the ride')
    else:
        run_ID = f'{args.dataset}_{args.model}'  
        print(f'Running without test mode active: THIS IS NOT A DRILL. \n Run ID: {run_ID}')

    paths = gnomadIC.setup_paths(run_ID)
    
    # Run chosen tasks
    gnomadIC.run_tasks(
        args.tasks, 
        paths = paths, 
        dataset = args.dataset,
        model = args.model,
        annotations = args.annotations,
        test = args.test, 
        controls = args.controls
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--targets', help='Path to target gene list',action='store')
    parser.add_argument('--test', help='Run tests',action='store_true',default=False)
    parser.add_argument('--controls',help='Include control genes',action='store_true',default=False)
    parser.add_argument('--dataset', help='Which dataset to use (one of gnomad, non_neuro, non_cancer, controls)', default='gnomad')
    parser.add_argument('--model', nargs= '+', help='Which model to apply (one of "standard", "syn_canonical", or "worst_csq" for now) - warning not implemented', default='standard')
    parser.add_argument('--annotations',help='Which annotations to apply (path to file)',action='store')
    parser.add_argument('--tasks', nargs='+', help='Which tasks to perform (select from download, model, summarise)')
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('-q','--quiet',help='Run in quiet mode',action='store_true',default=False)
    args = parser.parse_args()
    main(args)