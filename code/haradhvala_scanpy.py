import scanpy as sc


#sc_obj = sc.read_h5ad("/Users/dhasanova/Documents/ETH/HS23/data/haradhvala_scanpy/CART_fulldataset_clustered.h5ad")
#sc_obj.obs.to_csv('/Users/dhasanova/Documents/ETH/HS23/data/haradhvala_scanpy/md_scanpy.csv',sep=',')
#print(sc_obj.obs)

sc_obj = sc.read_h5ad("/Users/dhasanova/Documents/ETH/HS23/data/stator_results/CD8_baseline/output/unbinarised_cell_data.h5ad")
sc_obj.obs.to_csv('/Users/dhasanova/Documents/ETH/HS23/data/stator_results/CD8_baseline/Shiny/md_scanpy.csv',sep=',')
print(sc_obj.obs)