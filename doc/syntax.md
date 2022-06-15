# Adage Syntax

Here is a reference guide for the Adage functions. These are implemented in the `DAG` base class and the `Adage` derived class

## Root Operators

### `template<class O> void Initial(O op)`

Inbound root lambda. Allowed `O` lambdas are of the form:
- `std::function<void(Node*,NodePath*)>` (implement in your code as `auto op = [/*captures*/](Node *N, NodePath *bins){ /*body*/ }`, with no return value since it is `void`)
- `std::function<void(NodePath*,Node*)>`
- `std::function<void(NodePath*)>`
- `std::function<void(Node*)>`
- `std::function<void()>` (no arguments)

We will refer to this set of lambdas as the "standard" lambdas. Use the `Node` pointer to access the current control node in the traversal, and the `NodePath` pointer to access the current DAG path of nodes to this control node. See the respective header files and tutorials for how to use these objects in practice.

### `template<class O> void Final(O op)`

Outbound root lambda. Allowed `O` lambdas are standard (as listed above).


## Leaf Operators

### `template<class O> void Payload(O op)`

Payload operator, staged on the leaf node. Since there is no difference between inbound and outbound lambdas on the leaf node, only one staging template function is provided. This `Payload` function is meant to be used on `Adage` to act on the stored `Histos` objects.

Allowed `O` lambdas are of the form:
- `std::function<void(Histos*)>` (implement in your code as `auto op = [/*captures*/](Histos *H){ /*body*/ }`; this is the most common use case)
- `std::function<void(NodePath*)>`
- `std::function<void(Histos*,NodePath*)>`
- `std::function<void(NodePath*,Histos*)>`
- `std::function<void()>`

Similar to the standard set, all of these return `void`. The first lambda includes a `Histos` pointer; this is the most common use case, where you want to act on the `Histos` object on every multi-dimensional bin. If you need the bin path, include the `NodePath` pointer (though you can also access binning information from the `Histos` object itself).

### `template<class O> void LeafOp(O op)`

Low-level function to stage a standard lambda `O` on the leaf operator. It is likely you will find `Payload` much more useful. `LeafOp` is wrapped by `Payload` via a standard lambda function which calls `Adage::GetPayloadData(NodePath*)` to access the `Histos` object in the current multi-dimensional bin specified by the `NodePath` pointer.


## Subloops

### `template<class O1, class O2> void Subloop(std::vector<TString> layers, O1 opBefore, O2 opAfter)`
### `template<class O> void Subloop(std::vector<TString> layers, O opBefore)`
### `void Subloop(std::vector<TString> layers)`

These functions will make a subloop by creating a control node to control the layers listed in `std::vector<TString> layers`. The layers are iterated through in the specified order. The first `Subloop` template function is the full function, and the other two are "overloads", so you do not have to specify all (or any) of the lambdas. The inbound lambda will be `opBefore` and the outbound lambda will be `opAfter`.

The allowed `O1,O2,O` lambdas are the standard set, since at the subloop control node you only need access to the current node and/or its node path.

### `template<class O> void BeforeSubloop(std::vector<TString> layers, O op)`
### `template<class O> void AfterSubloop(std::vector<TString> layers, O op)`

These functions allow for an alternative specification of subloops. Again, `O` is a standard lambda.

### Conditional Controls

It can be useful to add conditional controls to subloop control nodes, to allow or disallow the iteration through the subloop. See [Adage documentation](adage.md) for more information.


## Multi-Payloads

### `template<class O1, class O2, class O3> void MultiPayload(std::vector<TString> layers, O1 opPayload, O2 opBefore, O3 opAfter)`
### `template<class O1, class O2> void MultiPayload(std::vector<TString> layers, O1 opPayload, O2 opBefore)`
### `template<class O1> void MultiPayload(std::vector<TString> layers, O1 opPayload)`

These three (overloaded) functions allow for the creation of multiple subloops, each with their own payload. The inbound and outbound lambdas, `opBefore` and `opAfter` are standard lambdas, and the payload lambda `opPayload` must have the same type you would use with the `Payload` staging function.

### `template<class O1, class O2, class O3> void MultiLeafOp(std::vector<TString> layers, O1 opLeaf, O2 opBefore, O3 opAfter)`

This is the low-level function wrapped by `MultiPayload`; it is not recommended to use this function directly. Behind the scenes, it will build the requested subloop, and then will call `Payload(opLeaf)` in the inbound lambda (along with `opBefore`), which overwrites any staged payload.
