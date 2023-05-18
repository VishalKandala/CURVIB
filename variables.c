#include "variables.h"

void initlist(List *ilist) {
  ilist->head = PETSC_NULL;
}

void insertnode(List *ilist, PetscInt Node)
{
  
  node *current,*new1;
  //node *new;
  current = ilist->head;

  PetscBool Exist = PETSC_FALSE;
  while(current) {
    if (Node == current->Node) {
      Exist = PETSC_TRUE;
    }
    if (Exist) break;
    current = current->next;
  }
  if (!Exist) {
    PetscMalloc(sizeof(node), &new1);
    new1->next = ilist->head;
    new1->Node = Node;
    ilist->head = new1;
  }
}

void destroy(List *ilist)
{
  node *current;
  while (ilist->head) {
    current = ilist->head->next;
    PetscFree(ilist->head);
    ilist->head = current;
  }
}

void InitIBMList(IBMList *ilist) {
  ilist->head = PETSC_NULL;
}

void AddIBMNode(IBMList *ilist, IBMInfo ibm_intp)
{
  IBMListNode *new1;
  PetscNew(&new1);
  new1->next = ilist->head;
  new1->ibm_intp = ibm_intp;
  ilist->head = new1;
}

void DestroyIBMList(IBMList *ilist)
{
  IBMListNode *current;
  while (ilist->head) {
    current = ilist->head->next;
    PetscFree(ilist->head);
    ilist->head = current;
  }
}

