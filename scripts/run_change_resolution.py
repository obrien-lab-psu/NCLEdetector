"""
python scripts/run_change_resolution.py 
--outdir TestingGrounds/CoarseGraining/ 
--pdbfile /path/to/structure.pdb 
--nscal 2 
--domain_file /path/to/domain_def.dat 
--ID test_protein
"""

def main(argv=None):

    import sys, os
    import argparse
    import time

    start_time = time.time()

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--outdir", type=str, required=True, help="output directory for results")
    parser.add_argument("--pdbfile", type=str, required=True, help="Path to the all-atom PDB file")
    parser.add_argument("--nscal", type=int, default=2, help="Coarse graining scale factor")
    parser.add_argument("--domain_file", type=str, required=True, help="Path to the domain definition file")
    parser.add_argument("--ID", type=str, required=True, help="Tag for the coarse graining process")
    args = parser.parse_args(argv)
    print(args)
    outdir = args.outdir
    pdbfile = args.pdbfile
    nscal = args.nscal
    domain_file = args.domain_file
    ID = args.ID

    # Import here so `run_change_resolution -h` works without OpenMM installed.
    from NCLEdetector.change_resolution import CoarseGrain, BackMapping

    ## Coarse grain the all-atom structure
    CoarseGrainer = CoarseGrain(outdir=outdir,
                                 ID=ID,
                                 pdbfile=pdbfile,
                                 nscal=nscal, 
                                 domain_file=domain_file)
    print(CoarseGrainer)
   
    # run the coarse graining for the all-atom pdbfile
    CGfiles = CoarseGrainer.run()
    print(CGfiles)

    # parse the prm and top file to make a OpenMM compatible .xml force field file
    CoarseGrainer.parse_cg_prm(prmfile=CGfiles['prm'], topfile=CGfiles['top'])


    ## Backmap the coarse grained structure to all-atom
    backMapper = BackMapping(outdir=outdir)
    print(f'BackMapper: {backMapper}')

    # backmap (using original pdbfile as aa_pdb reference)
    backMapper.backmap(cg_pdb=CGfiles['cor'], aa_pdb=pdbfile, TAG=ID)

    print(f'NORMAL TERMINATION - {time.time() - start_time} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
