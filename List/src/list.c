#include <stdio.h>
#include <stdlib.h>
#include "memory.h"
#include "list.h"
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/******************************************************************************
* Create a list object, initialize it, return it
******************************************************************************/
List *ListNew()
{
    List *lp;
    lp = SafeCalloc(1,sizeof(List));
    return lp;
}
/******************************************************************************
* Push an object to the back of the list
******************************************************************************/
void ListPushBack(List *lp,void *obj)
{
    ListNode *tmp_node=SafeCalloc(1,sizeof(ListNode));
    tmp_node->data=obj;
    tmp_node->next=NULL;
    if(lp->tail)
    {
        tmp_node->prev=lp->tail;
        lp->tail->next=tmp_node;
    }
    else
    {
        tmp_node->prev=NULL;
        lp->head=tmp_node;
    }
    lp->tail=tmp_node;
    lp->size++;
}



/******************************************************************************
* Push an object to the front of the list
******************************************************************************/
void ListPushFront(List *lp,void *obj)
{
    ListNode *tmp_node=SafeCalloc(1,sizeof(ListNode));
    tmp_node->data=obj;
    tmp_node->prev=NULL;
    if(lp->head)
    {
        tmp_node->next=lp->head;
        lp->head->prev=tmp_node;
    }
    else
    {
        tmp_node->next=NULL;
        lp->tail=tmp_node;
    }
    lp->head=tmp_node;
    lp->size++;
}


/******************************************************************************
* Erase the tail element from the list
******************************************************************************/
void *ListPopBack(List *lp)
{
    void *data_ptr;
    ListNode *tmp_node;
    switch(lp->size)
    {
    case 0:
        break;
    case 1:
        data_ptr = lp->tail->data;
        free(lp->tail);
        lp->head=NULL;
        lp->tail=NULL;
        lp->size=0;
        break;
    default:
        data_ptr = lp->tail->data;
        tmp_node=lp->tail->prev;
        free(lp->tail);
        lp->tail=tmp_node;
        tmp_node->next=NULL;
        (lp->size)--;
        break;
    }
    return data_ptr;
}



/******************************************************************************
* Erase the head element from the list
******************************************************************************/
void *ListPopFront(List *lp)
{
    void *data_ptr;
    ListNode *tmp_node;
    switch(lp->size)
    {
    case 0:
        break;
    case 1:
        data_ptr = lp->head->data;
        free(lp->head);
        lp->head=NULL;
        lp->tail=NULL;
        lp->size=0;
        break;
    default:
        data_ptr = lp->head->data;
        tmp_node=lp->head->next;
        free(lp->head);
        lp->head=tmp_node;
        tmp_node->prev=NULL;
        (lp->size)--;
        break;
    }
    return data_ptr;
}

/******************************************************************************
* Erase current list content
******************************************************************************/
void ListFree(List *lp,void free_func(void *))
{
    ListNode *tmp_node=lp->head;
    if(tmp_node)
    {
        tmp_node=tmp_node->next;
        while(tmp_node)
        {
            if(tmp_node->prev->data) free_func(tmp_node->prev->data);
            if(tmp_node->prev) free(tmp_node->prev);
            tmp_node=tmp_node->next;
        }
        tmp_node = lp->tail;
        if(tmp_node->data) free_func(tmp_node->data);
        if(tmp_node) free(tmp_node);
    }
    lp->head=NULL;
    lp->tail=NULL;
    lp->size=0;
}



/******************************************************************************
* Get the n-th element
* indexing starts from 0
******************************************************************************/
void *ListGetData(List *lp, unsigned int n)
{
    ListNode *node_ptr=lp->head;
    unsigned int i0;
    if(n>lp->size-1)
        return 0;
    for(i0 = 0; i0 < n; i0++, node_ptr = node_ptr->next) ;
    return node_ptr->data;
}
