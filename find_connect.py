import numpy as np

### Tests/examples ############################
def tests():
  ''' Examples/tests of use. '''
  set_tests()
  random_example()

def set_tests(size=10):
  print("No connections test...")
  test = np.eye(size,dtype=bool)
  connected_sets,home_node = find_connected_sets(test)
  print("result:")
  print(connected_sets)
  assert len(connected_sets) == size, 'No connection test fails.'
  print("...pass.")

  print("All connected test...")
  test = np.ones((size,size),dtype=bool)
  connected_sets,home_node = find_connected_sets(test)
  print("result:")
  print(connected_sets)
  assert len(connected_sets) == 1, 'All connection test fails.'
  print("...pass.")

  print("Set example test...")
  test = np.array([
      [1,0,0,1],
      [0,1,1,0],
      [0,1,1,0],
      [1,0,0,1] 
    ],dtype=bool)
  connected_sets,home_node = find_connected_sets(test)
  print("result:")
  print(connected_sets)
  assert len(connected_sets) == 2, 'All connection test fails.'
  print("...pass.")

def random_example(size=10,prob_connect=0.25):
  #np.random.seed(1003)
  test = np.random.rand(size,size)
  test = (test + test.T)/2. # I know it's not uniform, get off my back!
  test = test - np.eye(test.shape[0]) < prob_connect
  print("Random test:")
  print(test)

  connected_sets,home_node = find_connected_sets(test)
  print("Results:")
  print(connected_sets)

  print("Rearranged:")
  rearrange = arrange_to_blocks(test)
  print(rearrange)

### Utils  ################################

def arrange_to_blocks(symmatrix,tol=1e-4):
  ''' Rearrange a symmetric matrix to seperate blocks.

  Note: this can probably be generalized to nonsymmetric by taking 
  max of matrix and its transpose for each off-diagonal element.
  The problem is that find_connected_sets needs to be a metric.

  Args: 
    symmatrix (array): to be rearranged, in-place.
    tol (float): tolerance from blocks to be independent.
  '''
  ordering = connected_order(symmatrix,tol)

  symmatrix = symmatrix[ordering]
  symmatrix = symmatrix[:,ordering]

  return symmatrix

def connected_order(symmatrix,tol=1e-4):
  ''' Ordering that rearranges a symmetric matrix to seperate blocks.

  Note: this can probably be generalized to nonsymmetric by taking 
  max of matrix and its transpose for each off-diagonal element.
  The problem is that find_connected_sets needs to be a metric.

  Args: 
    symmatrix (array): to be rearranged, in-place.
    tol (float): tolerance from blocks to be independent.
  '''
  connections = np.abs(symmatrix) > tol
  connected_sets,_ = find_connected_sets(connections)

  ordering = []
  for head_node in connected_sets:
    ordering += connected_sets[head_node]
  return ordering

### Work horse ################################
def find_connected_sets(connects):
  ''' Find sets of nodes in graph that are connected.

  Args:
    connects (array-like): Bool array marking connected nodes.
  Returns:
    connected_set (dict): sets indexed by lowest-indexed node in the set.
    home_node (dict): for each node, the lowest-indexed node its connected to.
  '''

  # What is the lowest node each node is connected to?
  home_node = dict(zip(range(connects.shape[0]),range(connects.shape[1])))
  # What does each lowest node contain?
  connected_sets = dict(zip(range(connects.shape[0]), [set([i]) for i in range(connects.shape[0])] ))

  for node in range(connects.shape[0]):
    _update_connections(connects,node,home_node,connected_sets)
  
  return connected_sets,home_node

def _update_connections(connects,node,home_node,connected_sets):
  ''' Check neighbors for lower-index connection.'''
  nbrs = np.nonzero(connects[node])[0]
  for nbr in nbrs:
    # If they are in the same group: no updates needed.
    if home_node[node] == home_node[nbr]:
      pass
    # If this node trumps nbr, need to update nbr.
    elif home_node[node] < home_node[nbr]:
      _update_connections(connects,nbr,home_node,connected_sets)
    # If not, need to update this connected group.
    elif home_node[node] > home_node[nbr]:
      old_head = home_node[node]
      # All in this group now point to lower node.
      for connect_node in connected_sets[old_head]:
        home_node[connect_node] = home_node[nbr]
      # Merge connected group, and remove group from records.
      connected_sets[home_node[nbr]].update(connected_sets[old_head])
      del connected_sets[old_head]

if __name__ == '__main__':
  tests()
