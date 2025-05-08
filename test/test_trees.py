from unittest import TestCase

from ete3 import Tree

from dark.trees import ete3root


class NewickTreeTester:
    """
    Mixin for testing Newick trees.
    """

    def makeTree(self) -> Tree:
        return Tree(self.newick, format=1)  # type: ignore[attr-defined]

    def check(self, tree, expected):
        self.assertEqual(  # type: ignore[attr-defined]
            expected, tree.write(format=1, format_root_node=True)
        )


class Test_Reroot_Rooted_ABC(TestCase, NewickTreeTester):
    r"""
    Test the ete3root function on a rooted tree with three nodes, as below.

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
        result = ete3root(self.makeTree(), ("A",))
        self.check(result, "(A:1,(B:1,C:1)BC:1)ABC:0;")

    def testRoot_B(self) -> None:
        """
        If we root on the branch leading to B we should get the expected result.
        Note that the node with name BC is preserved but it is the node whose
        children are A/C after the re-rooting.
        """
        result = ete3root(self.makeTree(), ("B",))
        self.check(result, "(B:0.5,(C:1,A:2)BC:0.5)ABC:0;")

    def testRoot_C(self) -> None:
        """
        If we root on the branch leading to C we should get the expected result.
        Note that the node with name BC is preserved but it is the node whose
        children are A/C after the re-rooting.
        """
        result = ete3root(self.makeTree(), ("C",))
        self.check(result, "(C:0.5,(B:1,A:2)BC:0.5)ABC:0;")

    def testRoot_BC(self) -> None:
        """
        If we root on the branch leading to B/C we should get the same tree
        back since this would add a node on that branch and discard the old
        root, resulting in the same topology.
        """
        result = ete3root(self.makeTree(), ("B", "C"))
        self.check(result, "((B:1,C:1)BC:1,A:1)ABC:0;")

    def testDetach_A(self) -> None:
        """
        If we root on the branch leading to A and also ask for the outgroup
        to be detached, we should get back the B/C tree.
        """
        result = ete3root(self.makeTree(), ("A",), detach=True)
        self.check(result, "(B:1,C:1)BC:0;")

    def testDetach_BC(self) -> None:
        """
        If we root on the branch leading to B/C and also ask for the outgroup
        to be detached, we should get back the A tree.
        """
        result = ete3root(self.makeTree(), ("B", "C"), detach=True)
        self.check(result, "A:0;")

    def testRoot_C_and_scale(self) -> None:
        """
        If we root on the branch leading to C we should get the expected result,
        including when we scale the distances to the children of the new root.
        Note that the node with name BC is preserved but it is the node whose
        children are A/C after the re-rooting.
        """
        result = ete3root(self.makeTree(), ("C",), scale=0.25)
        # The distances to C and to the A/B node will both be 2 because the
        # original distance (of 1) to the C node is split between the two
        # children of the inserted root node (so 0.5 in each direction) and
        # we pass scale=0.25 which the new distances are divided by (and 0.5
        # / 0.25 = 2).
        self.check(result, "(C:2,(B:1,A:2)BC:2)ABC:0;")


class Test_Reroot_Unrooted_ABC(TestCase, NewickTreeTester):
    r"""
    Test the ete3root function on an unrooted tree with three nodes, as below.

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
        result = ete3root(self.makeTree(), ("A",))
        self.check(result, "(A:0.5,(B:1,C:1):0.5)ABC:0;")

    def testRoot_AB(self) -> None:
        """
        If we root on the branch leading to the MRCA of A/B we should get
        a tree with C on one side and an A/B subtree on the other. Note that
        the node with name ABC is preserved but it is the node whose children
        are B/C after the re-rooting.
        """
        result = ete3root(self.makeTree(), ("A", "B"))
        self.check(result, "(C:0.5,(A:1,B:1):0.5)ABC:0;")


class Test_Reroot_Unrooted_bigger_tree(TestCase, NewickTreeTester):
    r"""
    Test the ete3root function on the bigger unrooted tree below.

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
        result = ete3root(self.makeTree(), ("A",))
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
        self.check(result, "(A:0.5,((B:1,(C:1,D:1):1):1,(E:1,(F:1,G:1):1):1):0.5):0;")

    def testRoot_DG(self) -> None:
        """
        Re-rooting on the MRCA of D and G is the same as rooting on the branch that
        leads to A.
        """
        result = ete3root(self.makeTree(), ("D", "G"))
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
        self.check(result, "(A:0.5,((B:1,(C:1,D:1):1):1,(E:1,(F:1,G:1):1):1):0.5):0;")

    def testRoot_EG(self) -> None:
        """
        Re-rooting on the MRCA of E and G should work as expected.
        """
        result = ete3root(self.makeTree(), ("E", "G"))
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
        self.check(result, "((E:1,(F:1,G:1):1):0.5,(A:1,(B:1,(C:1,D:1):1):1):0.5):0;")

    def testRoot_CD(self) -> None:
        """
        Re-rooting on the MRCA of C and D should work as expected.
        """
        result = ete3root(self.makeTree(), ("C", "D"))
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
        self.check(result, "((C:1,D:1):0.5,(B:1,(A:1,(E:1,(F:1,G:1):1):1):1):0.5):0;")

    def testRoot_CD_detach(self) -> None:
        """
        Re-rooting on the MRCA of C and D should work as expected when we also then
        detach the C/D subtree.
        """
        result = ete3root(self.makeTree(), ("C", "D"), detach=True)
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
        self.check(result, "(B:1,(A:1,(E:1,(F:1,G:1):1):1):1):0;")
