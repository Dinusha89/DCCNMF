# DCCNMF
Matlab implementation for "Deep Complementary and Consensus Non-negative Matrix Factorization for multi-view clustering"

Matlab implementation of the following paper,
If you find it useful, please consider citing our work.

@article{gunawardena2024dccnmf,
  title={DCCNMF: Deep Complementary and Consensus Non-negative Matrix Factorization for multi-view clustering},
  author={Gunawardena, Sohan and Luong, Khanh and Balasubramaniam, Thirunavukarasu and Nayak, Richi},
  journal={Knowledge-Based Systems},
  volume={285},
  pages={111330},
  year={2024},
  publisher={Elsevier}
}

Abstract: Deep non-negative matrix factorization-based methods have recently been explored in multi-view clustering due to their ability to deal with complex non-linear data. Although a few methods exist that learn both complementary and consensus information simultaneously by adding new regularizations or adjusting hyperparameters without providing a solid architecture, they do not fully exploit this information. This paper proposes a novel method called Deep Complementary and Consensus Non-negative Matrix Factorization (DCCNMF) that combines the strengths of Non-negative Matrix Factorization and Coupled Non-negative Matrix Factorization using a novel architecture to simultaneously learn both complementary and consensus information present in multi-view data. Two manifold regularization terms, namely complementary manifold and consensus manifold are introduced to preserve the view-specific and view-shared geometric structures during the dimensionality reduction from higher to lower order. Further, smoothness regularization is employed to achieve distinct low-order representations by maximizing the cosine similarity among data points with similar orientations and minimizing it among data points with dissimilar orientations. DCCNMF and benchmark methods are evaluated on several real-world datasets and the results demonstrate that DCCNMF significantly outperforms state-of-the-art multi-view learning approaches.

# Demo
Run "Main.m" to see the provided example of BBCSport2views dataset.
