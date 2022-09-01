#ifndef _LIST_H
#define _LIST_H
//
// Encapsulation macros
#define ListNodeData(n,datatype) ((datatype)n->data)
#define ListNodeNext(n) (n->next)
#define ListNdoePrevious(n) (n->prev)
#define ListSize(l) (l->size)
#define ListHead(l) (l->head)
#define ListTail(l) (l->tail)
//
// Atomic type for list structure
// List will be run back and forth
typedef struct listnode
{
    void        *data;
    struct listnode *next;
    struct listnode *prev;
} ListNode;
// the list data structure
typedef struct list
{
    unsigned int    size;
    struct listnode *head;
    struct listnode *tail;
} List;
//
//
List *ListNew();
void ListPushFront(List *lp,void *obj);
void *ListPopFront(List *lp);
void ListPushBack(List *lp,void *obj);
void *ListPopBack(List *lp);
void *ListGetData(List *lp, unsigned int n);
void ListFree(List *lp,void free_func(void *));
//
#endif
