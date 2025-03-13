from WholeCellConsortiumModel import WholeCellConsortiumModel

wccm = WholeCellConsortiumModel()
wccm.generate_whole_network('test')
test_graph_2 = wccm.seek_optimal_network('test',['C00048'],['C00469','C00011'],2,2)
wccm.visualize_graph('test')