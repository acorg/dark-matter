from unittest import TestCase

from ete4 import Tree

from dark.trees import ete4root


class NewickTreeTester:
    """
    Mixin for testing Newick trees.
    """

    def makeTree(self) -> Tree:
        return Tree(self.newick, parser=1)  # type: ignore[attr-defined]

    def check(self, tree, expected):
        self.assertEqual(  # type: ignore[attr-defined]
            expected, tree.write(parser=1, format_root_node=True)
        )


class Test_Ancestor_MRCA_Unrooted(TestCase):
    """
    Some simple ancestor and MRCA tests on an unrooted tree.
    """

    newick = "(A, B, C);"

    def testSimple(self) -> None:
        tree = Tree(self.newick)
        a_nodes = list(tree.search_nodes(name="A"))
        assert len(a_nodes) == 1
        a_tip = a_nodes[0]
        ancestors = list(a_tip.ancestors(include_root=False))
        assert len(ancestors) == 1
        assert ancestors[0] == tree.root

    def testSimple_Parser1(self) -> None:
        tree = Tree(self.newick, parser=1)
        a_nodes = list(tree.search_nodes(name="A"))
        assert len(a_nodes) == 1
        a_tip = a_nodes[0]
        ancestors = list(a_tip.ancestors(include_root=False))
        assert len(ancestors) == 1
        assert ancestors[0] == tree.root

    def testAncestorsAndMRCAOfTwo(self) -> None:
        tree = Tree(self.newick)
        a_nodes = list(tree.search_nodes(name="A"))
        assert len(a_nodes) == 1
        a_tip = a_nodes[0]

        a_ancestors = list(a_tip.ancestors(include_root=False))
        assert len(a_ancestors) == 1
        assert a_ancestors[0] == tree.root

        b_nodes = list(tree.search_nodes(name="B"))
        assert len(b_nodes) == 1
        b_tip = b_nodes[0]

        b_ancestors = list(b_tip.ancestors(include_root=False))
        assert len(b_ancestors) == 1
        assert b_ancestors[0] == tree.root

        mrca = tree.common_ancestor(a_tip, b_tip)
        assert len(mrca.children) == 3
        assert mrca == tree.root

        included = {a_tip, b_tip}
        for tip in a_tip, b_tip:
            # The include_root is not important. Works either way (confirmed).
            ancestors = tip.ancestors(include_root=False)
            included.update(ancestors)

        excluded = set(mrca.children) - included
        assert len(excluded) == 1
        assert excluded.pop().name == "C"


class Test_Reroot_Rooted_ABC(TestCase, NewickTreeTester):
    r"""
    Test the ete4root function on a rooted tree with three nodes, as below.

       /-A
    --|
      |   /-B
       \-|
          \-C

    """

    newick = "(A, (B, C)BC)ABC;"

    def testRoot_A(self) -> None:
        """
        If we root on the branch leading to A we should get the same tree
        back since this would add a node on that branch and discard the old
        root, resulting in the same topology.
        """
        result = ete4root(self.makeTree(), ("A",))
        self.check(result, "(A,(B,C)BC);")

    def testRoot_B(self) -> None:
        """
        If we root on the branch leading to B we should get the expected result.
        Note that the node with name BC is preserved but it is the node whose
        children are A/C after the re-rooting.
        """
        result = ete4root(self.makeTree(), ("B",))
        self.check(result, "(B,(C,A)BC);")

    def testRoot_C(self) -> None:
        """
        If we root on the branch leading to C we should get the expected result.
        Note that the node with name BC is preserved but it is the node whose
        children are A/C after the re-rooting.
        """
        result = ete4root(self.makeTree(), ("C",))
        self.check(result, "(C,(B,A)BC);")

    def testRoot_BC(self) -> None:
        """
        If we root on the branch leading to B/C we should get the same tree
        back since this would add a node on that branch and discard the old
        root, resulting in the same topology.
        """
        result = ete4root(self.makeTree(), ("B", "C"))
        self.check(result, "((B,C)BC,A);")

    def testDetach_A(self) -> None:
        """
        If we root on the branch leading to A and also ask for the outgroup
        to be detached, we should get back the B/C tree.
        """
        result = ete4root(self.makeTree(), ("A",), detach=True)
        self.check(result, "(B,C)BC:0;")

    def testDetach_BC(self) -> None:
        """
        If we root on the branch leading to B/C and also ask for the outgroup
        to be detached, we should get back the A tree.
        """
        print(self.makeTree())
        result = ete4root(self.makeTree(), ("B", "C"), detach=True)
        self.check(result, "A:0;")


class Test_Reroot_Rooted_ABC_And_Scale(TestCase, NewickTreeTester):
    r"""
    Test the ete4root function (with subsequent scaling) on a rooted tree with
    three nodes, as below.

       /-A
    --|
      |   /-B
       \-|
          \-C

    """

    newick = "(A, (B, C:1)BC)ABC;"

    def testRoot_C(self) -> None:
        """
        If we root on the branch leading to C we should get the expected result,
        including when we scale the distances to the children of the new root.
        Note that the node with name BC is preserved but it is the node whose
        children are A/C after the re-rooting.
        """
        result = ete4root(self.makeTree(), ("C",), scale=0.25)
        # The distances to C and to the A/B node will both be 2 because the
        # original distance (of 1) to the C node is split between the two
        # children of the inserted root node (so 0.5 in each direction) and
        # we pass scale=0.25 which the new distances are divided by (and 0.5
        # / 0.25 = 2).
        # self.check(result, "(C:2,(B:1,A:2)BC:2)ABC:0;")
        self.check(result, "(C:2,(B,A)BC:2);")

        result = ete4root(self.makeTree(), ("C",), scale=0.1)
        self.check(result, "(C:5,(B,A)BC:5);")

        result = ete4root(self.makeTree(), ("C",), scale=2.0)
        self.check(result, "(C:0.25,(B,A)BC:0.25);")


class Test_Reroot_Unrooted_ABC(TestCase, NewickTreeTester):
    r"""
    Test the ete4root function on an unrooted tree with three nodes, as below.

       /-A
      |
    --|--B
      |
       \-C

    """

    newick = "(A, B, C)ABC;"

    def testRoot_A(self) -> None:
        """
        If we root on the branch leading to A we should get the expected result.
        Note that the node with name ABC is preserved but it is the node whose
        children are B/C after the re-rooting.
        """
        result = ete4root(self.makeTree(), ("A",))
        self.check(result, "(A,(B,C)ABC);")

    def testRoot_AB(self) -> None:
        """
        If we root on the branch leading to the MRCA of A/B we should get
        a tree with C on one side and an A/B subtree on the other. Note that
        the node with name ABC is preserved but it is the node whose children
        are B/C after the re-rooting.
        """
        result = ete4root(self.makeTree(), ("A", "B"))
        self.check(result, "(C,(A,B)ABC);")


class Test_Reroot_Unrooted_bigger_tree(TestCase, NewickTreeTester):
    r"""
    Test the ete4root function on the bigger unrooted tree below.

       /-A
      |
      |   /-B
      |--|
    --|  |   /-C
      |   \-|
      |      \-D
      |
      |   /-E
       \-|
         |   /-F
          \-|
             \-G

    """

    newick = "(A, (B, (C, D)), (E, (F, G)));"

    def testRoot_A(self) -> None:
        """
        Re-rooting on the branch that leads to A should work as expected.
        """
        result = ete4root(self.makeTree(), ("A",))
        # result looks like this (you can see it via print(result)).
        r"""
           /-A
          |
          |      /-B
        --|   /-|
          |  |  |   /-C
          |  |   \-|
           \-|      \-D
             |
             |   /-E
              \-|
                |   /-F
                 \-|
                    \-G
        """
        self.check(result, "(A,((B,(C,D)),(E,(F,G))));")

    def testRoot_DG(self) -> None:
        """
        Re-rooting on the MRCA of D and G is the same as rooting on the branch that
        leads to A.
        """
        result = ete4root(self.makeTree(), ("D", "G"))
        # result looks like this (you can see it via print(result)).
        r"""
           /-A
          |
          |      /-B
        --|   /-|
          |  |  |   /-C
          |  |   \-|
           \-|      \-D
             |
             |   /-E
              \-|
                |   /-F
                 \-|
                    \-G
        """
        self.check(result, "(A,((B,(C,D)),(E,(F,G))));")

    def testRoot_EG(self) -> None:
        """
        Re-rooting on the MRCA of E and G should work as expected.
        """
        result = ete4root(self.makeTree(), ("E", "G"))
        # result looks like this (you can see it via print(result)).
        r"""
              /-E
           /-|
          |  |   /-F
          |   \-|
        --|      \-G
          |
          |   /-A
           \-|
             |   /-B
              \-|
                |   /-C
                 \-|
                    \-D
        """
        self.check(result, "((E,(F,G)),(A,(B,(C,D))));")

    def testRoot_CD(self) -> None:
        """
        Re-rooting on the MRCA of C and D should work as expected.
        """
        result = ete4root(self.makeTree(), ("C", "D"))
        # result looks like this (you can see it via print(result)).
        r"""
              /-C
           /-|
          |   \-D
        --|
          |   /-B
           \-|
             |   /-A
              \-|
                |   /-E
                 \-|
                   |   /-F
                    \-|
                       \-G
        """
        self.check(result, "((C,D),(B,(A,(E,(F,G)))));")

    def testRoot_CD_detach(self) -> None:
        """
        Re-rooting on the MRCA of C and D should work as expected when we also then
        detach the C/D subtree.
        """
        result = ete4root(self.makeTree(), ("C", "D"), detach=True)
        # result looks like this (you can see it via print(result)).
        r"""
           /-B
        --|
          |   /-A
           \-|
             |   /-E
              \-|
                |   /-F
                 \-|
                    \-G
        """
        self.check(result, "(B,(A,(E,(F,G)))):0;")
