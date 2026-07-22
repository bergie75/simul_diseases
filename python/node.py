import numpy as np
from copy import deepcopy

class Node:
    def __init__(self, index, transmit_prob=[1,1], fall_ill_prob=[1,1], recover_prob=[1,1],
                  status=["S", "S"], neighbors=[], blocking_disease=None,
                  inf_decision_weights=[1,1], prior_resistance=[0,0]):
        # random number generator to handle all node instances
        self.rng = np.random.default_rng()

        self.index = index  # holds index metadata, works with a list of all nodes
        self.transmit_prob = np.array(transmit_prob)  # probability of transmitting each disease to 
        self.fall_ill_prob = np.array(fall_ill_prob)  # probability of moving from exposed to infected
        self.recover_prob = np.array(recover_prob)  # probability of moving from infected to recovered
        self.status = deepcopy(status)  # current compartment out of SEIR for each disease
        self.neighbors = neighbors  # list of all neighboring nodes
        self.blocking_disease = blocking_disease  # prevents further exposures if set to a disease index
        self.inf_decision_weights = inf_decision_weights  # how easily one disease can prevent the other from getting a foothold
        self.prior_resistance = prior_resistance

        # temp variables
        self.temp_exposure_queue = np.array([0.0,0.0])
    
    def update_neighbors(self, new_neighbors):
        self.neighbors = new_neighbors
    
    def accept_infection(self, disease_index, day=None):
        status_readout = [False, disease_index, self.index, self.status, day, False]
        if self.blocking_disease is None and self.status[disease_index]=="S" and self.rng.uniform()>self.prior_resistance[disease_index]:
            self.temp_exposure_queue[disease_index] += self.inf_decision_weights[disease_index]
        elif self.blocking_disease is not None:
            status_readout[0] = True  # an exposure has been blocked
            if self.blocking_disease == disease_index:
                status_readout[-1] = True  # this is a "self-block"
        
        return status_readout
    
    def choose_exposure_event(self):
        # check blocking disease first for efficiency
        if self.blocking_disease is None:
            total_disease_weight = np.sum(self.temp_exposure_queue)
            if total_disease_weight > 0:
                dis_zero_prob = self.temp_exposure_queue[0]/total_disease_weight
                chosen_disease = self.rng.binomial(1, 1-dis_zero_prob)
                self.temp_exposure_queue *= 0.0  # resets exposure queue
                self.blocking_disease = chosen_disease
                self.status[chosen_disease] = "E"
                return chosen_disease
    
    def progress_status(self):
        for disease_index, current_status in enumerate(self.status):
            if current_status == "E" and self.rng.uniform()<= self.fall_ill_prob[disease_index]:
                self.status[disease_index] = "I"
            elif current_status == "I" and self.rng.uniform() <= self.recover_prob[disease_index]:
                self.status[disease_index] = "R"
                self.blocking_disease = None  # disease releases hold on the node, allowing new diseases to infect
    
    def transmit(self, disease_index, day=None, trans_modifier=1.0):
        # check if currently exposed to the disease. If yes, transmit probabilistically to neighbors
        # if not, pass
        blocking_events = []
        if self.status[disease_index] == "E":
            for neighbor in self.neighbors:
                if self.rng.uniform() <= self.transmit_prob[disease_index]*trans_modifier:
                    blocking_events.append(neighbor.accept_infection(disease_index, day=day))
        
        return blocking_events

if __name__ == "__main__":
    # test node behaviors
    recip1 = Node(2)
    recip2 = Node(3)
    print("initial status check")
    print(recip1.status)
    print(recip2.status)
    totally_contagious = Node(1, neighbors=[recip1, recip2], status=["E", "E"])
    totally_contagious.transmit(0)
    print("prior to accept transmit")
    print(recip1.status)
    print(recip2.status)
    recip1.choose_exposure_event()
    recip2.choose_exposure_event()
    print("post acceptance")
    print(recip1.status)
    print(recip2.status)
    print("progress status 1")
    recip1.progress_status()
    print(recip1.status)
    print(recip1.blocking_disease)
    print("attempt transmit second disease")
    totally_contagious.transmit(1)
    recip1.choose_exposure_event()
    recip2.choose_exposure_event()
    print(recip1.status)
    print(recip1.blocking_disease)
    print("post progress 2")
    recip1.progress_status()
    print(recip1.status)
    print(recip1.blocking_disease)