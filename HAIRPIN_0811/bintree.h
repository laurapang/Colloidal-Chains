typedef struct Node {
  double key;
  int bead1;
  int bead2;
  struct Node* left;
  struct Node* right;
} Node;

typedef struct Tree {
  unsigned short size;
  Node* head;
} Tree;

void initTree( Tree* t );

void initNode( Node* n, double key, int bead1, int bead2 );

void add( Node* current, Node* n );

void addNode( Tree* t, double key, int bead1, int bead2 );

Node* addInorder( Node* n, Node* list );

Node* getInorder( Tree* t );

void destroyNode( Node* n );

void destroyTree( Tree* t );

// Test functions

void inorder( Node* n );

void printInorder( Tree* t );