#include <stdio.h>
#include <stdlib.h>
#include "bintree.h"

// Useful functions

// Tree constructor
void initTree( Tree* t ) {
  t->size = 0;
  t->head = NULL;
}

// Node constructor
void initNode( Node* n, double key, int bead1, int bead2 ) {
  n->key = key;
  n->bead1 = bead1;
  n->bead2 = bead2;
  n->left = NULL;
  n->right = NULL;
}

// Recursively add a node to another or one of its children
void add( Node* current, Node* n ) {
  if ( n->key <= current->key ) {
    if ( !current->left ) // if NULL
      current->left = n;
    else
      add( current->left, n );
  } else {
    if ( !current->right )
      current->right = n;
    else
      add( current->right, n );
  }
}

// Make a new node and add it to the tree
void addNode( Tree* t, double key, int bead1, int bead2 ) {
  Node* n = ( Node* ) malloc( sizeof( Node ) );
  initNode( n, key, bead1, bead2 );

  if ( !t->head )
    t->head = n;
  else
    add( t->head, n );

  t->size++;
}

// Helper function for recursively building inorder list
// Every time a node is added the pointer in the list
// is incremented, pointing to the next element
// in the list
Node* addInorder( Node* n, Node* list ) {
  if ( n->left ) // if not NULL
    list = addInorder( n->left, list );
  *list = *n;
  list++;
  if ( n->right )
    list = addInorder( n->right, list );
  return list;
}

// Get array of nodes in ascending order by key
Node* getInorder( Tree* t ) {
  Node* nodes = ( Node* ) malloc( t->size * sizeof( Node ) );

  // Get ptr to front of list so we dont lose the array
  Node* ptr = nodes;
  (void) addInorder( t->head, ptr );
  
  return nodes;
}

void destroyNode( Node* n ) {
  if ( n->left != NULL )
    destroyNode( n->left );
  if ( n->right != NULL )
    destroyNode( n->right );
  free( n );
}

void destroyTree( Tree* t ) {
  destroyNode( t->head );
  free( t );
}

// Test functions

// Prints left child, parent, then right child recursively
void inorder( Node* n ) {
  if ( n->left )
    inorder( n->left );
  printf( "%.4f\t", n->key );
  if ( n->right )
    inorder( n->right );
}

// Prints keys of tree in ascending order on a single line
void printInorder( Tree* t ) {
  inorder( t->head );
  printf( "\n" );
}

// Example of using the binary tree
#if defined(TESTING)

int main() {

  // Init the tree and fill it with nodes with arbitrary keys
  Tree* t = (Tree*) malloc( sizeof( Tree ) );
  initTree( t );

  double numbers[7] = { 43.0, 21.0, 12.0, 100.0, 72.0, 81.0, 101.0 };

  for ( unsigned char i = 0; i < 7; i++ )
    addNode( t, numbers[i], 0, 0);

  // Test that the tree initialization worked
  printInorder(t);

  // Retrieve an array of the nodes in order and print their keys
  // an array (a pointer) of node structs
  Node* inordered = getInorder( t );

  for ( unsigned char i = 0; i < t->size; i++ )
    printf( "%.4f\t", inordered[i].key );    
  // pointers to structs use '->', normal structs use '.' like java objects
  printf( "\n" );

  // Free the dynamically allocated memory
  // and make sure pointers dont point there anymore
  destroyTree( t );
  t = NULL;
  free( inordered );
  inordered = NULL;

  return 0;
}

#endif
