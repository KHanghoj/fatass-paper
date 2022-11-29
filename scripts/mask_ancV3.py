import sys
import argparse
import numpy as np

def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "--loglike",
        type=str,
        required=True,
        help=""
    )

    parser.add_argument(
        "--decoding",
        type=str,
        required=False,
        help=""
    )

    parser.add_argument(
        "--posterior",
        type=str,
        required=False,
        help=""
    )

    parser.add_argument(
        "--outbase",
        type=str,
        required=True,
        help=""
    )

    parser.add_argument(
        '--het_mask',
        action='store_true',
        help="Masking per haplotype"
    )
    
    parser.add_argument(
        '--max_missing',
        default=0.95,
        type=float,
        help="Samples with more missingness will be excluded. [0,1] Set to 1 to disable."
    )

    parser.add_argument(
        '--posterior_cutoff',
        default=0.95,
        type=float,
        help="Samples with posterior les than will be masked."
    )

    parser.add_argument(
        '--conservative',
        action='store_true',
        help="set windows to missing if next to a change in ancestry (development)"
    )

    parser.add_argument(
        "--labels",
        required=True,
        type=str,
    )

    parser.add_argument(
        "--names",
        required=True,
        type=str,
    )



    return parser.parse_args(args=None if argv else ['--help'])


def main(argv):

    args = parse_args(argv)
    if args.conservative:
        print("--convservative not implemented", file=sys.stderr)
        sys.exit(1)
    
    if (args.posterior and args.decoding) or ((not args.posterior) and (not args.decoding)):
        print("Provide either --posterior or --decoding. Not both and not none of them :)", file=sys.stderr)
        sys.exit(1)

    log = open(f"{args.outbase}.log", 'w')

    for arg, value in vars(args).items():
        print(f"Argument {arg}={value}", file=log)

    if args.labels:
        print("parsing labels", file=log)
        labels = list(np.loadtxt(args.labels, dtype=str))
    if args.names:
        print("parsing names", file=log)
        names = list(np.loadtxt(args.names, dtype=str))

    if args.het_mask:
        print("masking per haplotype", file=log)
    else:
        print("masking homozygous ancestry", file=log)
        
    if args.posterior:
        post = np.exp(np.load(args.posterior))
        K = post.shape[2]
    elif args.decoding:
        ## loading decoding
        post = np.loadtxt(args.decoding, dtype=int).T
        K = post.max()+1
    else:
        print("should not happen", file=fh)
        sys.exit(1)
    # path = np.loadtxt(args.decoding, skiprows=1, dtype=int).reshape(Ws, N, HAPLOTYPES)
    Ls = [np.load(args.loglike) for _ in range(K)]
    Ws,Hs,Cs = Ls[0].shape
    print(Ls[0].shape)
    print(len(names))
    print(len(labels))
    N = Hs//2
    HAPLOTYPES=2
    assert N == len(names)
    assert N == len(labels)

    
    names_k = []
    labels_k = []
    for k in range(K):
        names_k.extend([n+f"K{k}" for n in names])
        labels_k.extend([n+f"K{k}" for n in labels])

    missing = np.log(np.ones(Cs)/Cs)
    window_idx1 = list(range(Ws-1))
    window_idx2 = list(range(1,Ws))

    for x in range(N):
        
        # Two haplotypes
        idx1 = x*2
        idx2 = x*2+1
        # if ind_unadmixed[x]:
        #     continue
        for k in range(K):
            if not args.het_mask:
                if args.posterior:
                    mask1 = post[:, idx1, k] < args.posterior_cutoff
                    mask2 = post[:, idx2, k] < args.posterior_cutoff
                else:
                    mask1 = post[:, idx1]!=k
                    mask2 = post[:, idx2]!=k
                mask_clusters1 =  np.logical_or(mask1, mask2)
                Ls[k][mask_clusters1, idx1, :] = missing
                Ls[k][mask_clusters1, idx2, :] = missing
            else:
                for idx in [idx1, idx2]:
                    if args.posterior:
                        mask = post[:, idx, k] < args.posterior_cutoff
                    else:
                        mask = post[:, idx] != k
                    Ls[k][mask, idx, :] = missing
                

        #     if args.conservative:

        #         ## haplotype 1
        #         exclude_temp = np.where(np.not_equal(path[window_idx1,x,0], path[window_idx2,x,0]))[0]
        #         exclude1 = np.concatenate((exclude_temp, exclude_temp+1))

        #         ## haplotype 2
        #         exclude_temp = np.where(np.not_equal(path[window_idx1,x,1], path[window_idx2,x,1]))[0]
        #         exclude2 = np.concatenate((exclude_temp, exclude_temp+1), axis=None)

        #         ## removes transistions if observed on either haplotypes
        #         exclude = np.concatenate((exclude1, exclude2))
        #         exclude = np.unique(exclude)

        #         L[exclude, idx1, :] = missing
        #         L[exclude, idx2, :] = missing
        #         L2[exclude, idx1, :] = missing
        #         L2[exclude, idx2, :] = missing

        # else:
        #     ## one ancestry
        #     mask_clusters1 =  mask1[:,0]<args.posterior_cutoff
        #     L[mask_clusters1, idx1, :] = missing

        #     mask_clusters2 =  mask2[:,0]<args.posterior_cutoff
        #     L[mask_clusters2, idx2, :] = missing

        #     ## other ancestry
        #     mask_clusters1 =  mask1[:,1]<args.posterior_cutoff
        #     L2[mask_clusters1, idx1, :] = missing

        #     mask_clusters2 =  mask2[:,1]<args.posterior_cutoff
        #     L2[mask_clusters2, idx2, :] = missing

        #     ## conservative
        #     if args.conservative:
        #         exclude_temp = np.where(np.not_equal(path[window_idx1,x,0], path[window_idx2,x,0]))[0]
        #         exclude1 = np.concatenate((exclude_temp, exclude_temp+1))
        #         L[exclude1, idx1, :] = missing
        #         L2[exclude1, idx1, :] = missing

        #         exclude_temp = np.where(np.not_equal(path[window_idx1,x,1], path[window_idx2,x,1]))[0]
        #         exclude2 = np.concatenate((exclude_temp, exclude_temp+1))
        #         L[exclude2, idx2, :] = missing
        #         L2[exclude2, idx2, :] = missing

    test = np.concatenate(Ls, axis=1)
    print(test.shape)
    newHs = test.shape[1]
    abs2 = np.all(np.isclose(test, np.log(1.0/Cs)), axis=2).reshape(Ws, newHs//2, 2)
    missing_prop = np.mean(abs2, axis=(0,2))
    with open(f"{args.outbase}.missingness", 'w') as fh:
        for n,l,m in zip(names_k, labels_k, missing_prop):
            print(n,l,m, file=fh)
    b = missing_prop<=args.max_missing
    keep = np.repeat(b,2)
    np.save(f"{args.outbase}.npy", test[:,keep,:])

    with open(f"{args.outbase}.labels", 'w') as fh:
        for idx, label in enumerate(labels_k):
            if b[idx]:
                print(label, file=fh)

    with open(f"{args.outbase}.names", 'w') as fh:
        for idx, name in enumerate(names_k):
            if b[idx]:
                print(name, file=fh)


    log.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
