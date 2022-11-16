
class Node:

    __NodeIdx = 0

    def __init__(self, Data, ParentNode=None, Prior=0.0):
        Node.__NodeIdx       += 1
        self.__NodeIdx        = Node.__NodeIdx
        self.Data             = Data
        self.ParentNode       = ParentNode
        self.ChildNodes       = []
        self.NodeVisited      = 0
        self.ActionValue      = 0.0
        self.PriorProbability = Prior
        self.Reward           = 0.0

    def GetNodeIdx(self):
        return self.__NodeIdx

    def GetNodeData(self):
        return self.Data

    def GetParentNode(self):
        return self.ParentNode

    def GetChildNodes(self):
        return self.ChildNodes

    def GetNumVisited(self):
        return self.NodeVisited

    def GetActionValue(self):
        return self.ActionValue

    def GetPriorProbability(self):
        return self.PriorProbability

    def GetReward(self):
        return self.Reward

    def IncreaseNodeVisit(self):
        self.__NodeIdx += 1

    def AddChildNode(self, childnode):
        self.ChildNodes.append(childnode)

    def SetRewardScore(self, reward):
        self.reward = reward

    def SetActionValue(self, action):
        self.ActionValue = action

class Tree:

    def __init__(self, Data):
        self.RootNode = Node(Data)
        self.TreeNodes = [[self.RootNode]]

    def GetTreeNode(self, TreeNodeIdx):
        for i, x in enumerate(self.TreeNodes):
            for j, node in enumerate(x):
                if TreeNodeIdx == node.GetNodeIdx():
                    return i, j, node

    def AddChildNode(self, ParentNode, ChildNodeData, Prior=None):
        i, j, _ = self.GetTreeNode(ParentNode.GetNodeIdx())
        ChildNode = Node(ChildNodeData, ParentNode, Prior)
        if len(self.TreeNodes) == i + 1:
            self.TreeNodes.append([])
        self.TreeNodes[i+1].append(ChildNode)
        self.TreeNodes[i][j].AddChildNode(ChildNode)

    def IsRootNode(self, TreeNode):
        if TreeNode.GetNodeIdx() == self.RootNode.GetNodeIdx():
            return True
        else:
            return False

    def IsLeafNode(self, TreeNode):
        if not TreeNode.GetChildNodes():
            return True
        else:
            return False

    def GetRootNode(self):
        return self.RootNode

    def GetLeafNodes(self):
        return [node for x in self.TreeNodes for node in x if self.IsLeafNode(node)]

    def GetBranchs(self, ParentNode):
        Lv, _, _  = self.GetTreeNode(ParentNode.GetNodeIdx())
        LeafNodes = self.GetLeafNodes()
        Branchs   = []
        for node in LeafNodes:
            LeafLv, _, _ = self.GetTreeNode(node.GetNodeIdx())
            if LeafLv < Lv:
                continue
            branch = [node]
            while LeafLv > Lv:
                if not self.IsRootNode(node):
                    node = node.GetParentNode()
                    branch.append(node)
                    LeafLv -= 1
            if branch[-1].GetNodeIdx() == ParentNode.GetNodeIdx():
                Branchs.append(branch)
        return Branchs

class MCTS:

    def __init__(self, Target, EndState, SearchIter=10, SimIteration=100):
        self.TreeSearch = Tree(Target)
        self.EndState   = EndState
        self.SearchIter = SearchIter
        self.SimIter    = SimIteration

    def __Selection(self):
        pass

    def __Expansion(self):
        pass

    def Simulation(self):
        pass

    def __Update(self):
        pass

