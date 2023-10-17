from tqdm import tqdm
from .kmer import KMer
import matplotlib.pyplot as plt

# Choose kmer
def choose_k(kmax, problem):

    accumulated_metrics = {
        'train_accuracy' : [],
        'test_accuracy' : [],
        'b_accuracy' : [],
        'f1' : [],
        'auc' : []
    }

    for k in tqdm(range(1, 11)):
        kmer = KMer(k = k, problem = problem)
        X = kmer.make_embedding()
        X = kmer.normalize_features()
        X = kmer.scalar_embedding(100)

        metrics = kmer.k_fold(X = kmer.X, y = kmer.y, k = 'loocv', plot = False)

        for metric in metrics:
                    accumulated_metrics[metric].append(metrics[metric])

    # Plot
    for metric, values in accumulated_metrics.items():
        plt.plot(list(range(1, kmax + 1)), values, label = metric)

    plt.ylabel('Metric value')
    plt.xlabel('k')
    plt.legend()
    plt.title(kmer.title)
    plt.show()

def main():
    '''Program flow.'''

    # 1)
    #choose_k(10, problem = 'SS')
    
    # 2)
    kmer = KMer(k = 2, problem = 'C1C4')
    X = kmer.make_embedding()
    X = kmer.normalize_features()
    X = kmer.scalar_embedding(100)

    
    # 3)
    #metrics = kmer.range_kfold(X = kmer.X, y = kmer.y)
    #print(metrics)

    # 4)
    print(kmer.k_fold(X = kmer.X, y = kmer.y, k = 'loocv', plot = True))


    # 5)
    #important_features = kmer.feature_importance()
    #kmer.kfold_per_feature(k = 'loocv')
    #print(kmer.n_features_for_metric('b_accuracy', 1))

    
    # 6)
    #kmer.most_similar_cluster()

    
    # 7)
    #print('Important features')
    #important_features = kmer.feature_importance()
    #print('Map features')
    #mapped_features = kmer.map_features(33)
    #print('Generate report')
    #kmer.features_report(33)
    #print('To PyMol')
    #kmer.report2pymol(33)

    # 8)
    #important_features = kmer.feature_importance()
    #table = kmer.map_features_report(34)
    #table.to_csv('table2.csv')

if __name__ == '__main__':
    main()