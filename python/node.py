import numpy as np
from copy import deepcopy
import yaml
import os

cwd = os.getcwd()
par_folder = os.path.join(cwd, "parameters")

# locate configuration file
with open(os.path.join(par_folder, "config.yaml"), "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)

# unpack variables to use in default Node generator
def_trans_prob = cfg["trans_prob"]
def_fall_ill = cfg["fall_ill"]
def_dis_weight = cfg["dis_weight"]
def_prior_res = cfg["prior_res"]
def_large_resp = cfg["large_resp"]

# separate out, requires a calculation
rec_weights = cfg["rec_weights"]
num_days = cfg["num_days"]
def_rec_prob = [x/num_days for x in rec_weights]

class Node:
    def __init__(self, index, trans_prob=def_trans_prob, fall_ill=def_fall_ill, rec_prob=def_rec_prob,
                  status=["S", "S"], neighbors=[], blocking_disease=None,
                  dis_weight=def_dis_weight, prior_res=def_prior_res, large_resp=def_large_resp):
        # random number generator to handle all node instances
        self.rng = np.random.default_rng()

        self.index = index  # holds index metadata, works with a list of all nodes
        self.trans_prob = np.array(trans_prob)  # probability of transmitting each disease to 
        self.fall_ill = np.array(fall_ill)  # probability of moving from exposed to infected
        self.rec_prob = np.array(rec_prob)  # probability of moving from infected to recovered
        self.status = deepcopy(status)  # current compartment out of SEIR for each disease
        self.neighbors = neighbors  # list of all neighboring nodes
        self.blocking_disease = blocking_disease  # prevents further exposures if set to a disease index
        self.dis_weight = np.array(dis_weight)  # how easily one disease can prevent the other from getting a foothold
        self.prior_res = np.array(prior_res)
        self.large_resp = large_resp

        # temp variables
        self.temp_exposure_queue = np.array([0.0,0.0])
        self.caused_large_resp = [False, False]
    
    def set_block(self, blocking_disease):
        self.blocking_disease = blocking_disease

    def set_status(self, disease_index, disease_status):
        self.status[disease_index] = disease_status

    def reset_temps(self):
        self.temp_exposure_queue = np.array([0.0, 0.0])
        self.caused_large_resp = [False, False]

    def dict_update(self, config_opts):
        for var in config_opts.keys():
            setattr(self, var, deepcopy(config_opts[var]))

    def update_neighbors(self, new_neighbors):
        self.neighbors = new_neighbors
    
    def accept_infection(self, disease_index, day=None):
        status_readout = [False, disease_index, self.index, deepcopy(self.status), day, False]
        if self.blocking_disease is None and self.status[disease_index]=="S" and self.rng.uniform()>self.prior_res[disease_index]:
            self.temp_exposure_queue[disease_index] += self.dis_weight[disease_index]
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
                self.status[chosen_disease] = "I"
                self.caused_large_resp[chosen_disease] = (self.rng.uniform() <= self.large_resp[chosen_disease])
                return chosen_disease
    
    def progress_status(self):
        for disease_index, current_status in enumerate(self.status):
            if current_status == "I" and self.rng.uniform()<= self.fall_ill[disease_index]:
                if self.caused_large_resp[disease_index]:
                    self.status[disease_index] = "H"
                else:
                    self.status[disease_index] = "R"
                    self.blocking_disease = None
            elif current_status == "H" and self.rng.uniform() <= self.rec_prob[disease_index]:
                self.status[disease_index] = "R"
                self.blocking_disease = None  # disease releases hold on the node, allowing new diseases to infect
    
    def transmit(self, nodelist, disease_index, day=None, trans_mod=1.0):
        # check if currently exposed to the disease. If yes, transmit probabilistically to neighbors
        # if not, pass
        blocking_events = []
        if self.status[disease_index] == "I":
            for neighbor in self.neighbors:
                if self.rng.uniform() <= self.trans_prob[disease_index]*trans_mod:
                    blocking_events.append(nodelist[neighbor].accept_infection(disease_index, day=day))
        
        return blocking_events

if __name__ == "__main__":
    new_node = Node(12)
    new_node.dict_update({"status": ["I", "R"], "blocking_disease": 0})
    print(new_node.status)
    print(new_node.blocking_disease)
    # test node behaviors
    # recip1 = Node(2)
    # recip2 = Node(3)
    # print("initial status check")
    # print(recip1.status)
    # print(recip2.status)
    # totally_contagious = Node(1, neighbors=[recip1, recip2], status=["I", "I"])
    # totally_contagious.transmit(0)
    # print("prior to accept transmit")
    # print(recip1.status)
    # print(recip2.status)
    # recip1.choose_exposure_event()
    # recip2.choose_exposure_event()
    # print("post acceptance")
    # print(recip1.status)
    # print(recip2.status)
    # print("progress status 1")
    # recip1.progress_status()
    # print(recip1.status)
    # print(recip1.blocking_disease)
    # print("attempt transmit second disease")
    # totally_contagious.transmit(1)
    # recip1.choose_exposure_event()
    # recip2.choose_exposure_event()
    # print(recip1.status)
    # print(recip1.blocking_disease)
    # print("post progress 2")
    # recip1.progress_status()
    # print(recip1.status)
    # print(recip1.blocking_disease)