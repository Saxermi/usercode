import torch

def test_cuda():
    if torch.cuda.is_available():
        num_gpus = torch.cuda.device_count()
        info = [f"CUDA is available. {num_gpus} GPU(s) detected."]
        for i in range(num_gpus):
            name = torch.cuda.get_device_name(i)
            info.append(f"GPU {i}: {name}")
        return "\n".join(info)
    else:
        return "CUDA is NOT available."

if __name__ == "__main__":
    print(test_cuda())
